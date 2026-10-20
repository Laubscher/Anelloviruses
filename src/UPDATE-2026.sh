#!/usr/bin/env bash
# 2026: updates for the six small genera.
# Requires Python >= 3.8. The update action also requires MUSCLE 3 and ORF1.sh.

set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
exec "${PYTHON:-python3}" - "$SCRIPT_DIR" "${BASH_SOURCE[0]}" "$@" <<'PYTHON'

import argparse
import csv
import datetime as dt
import fcntl
from fractions import Fraction
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request

GENERA = {
    "he": "Hetorquevirus", "sade": "Sadetorquevirus",
    "yod": "Yodtorquevirus", "lamed": "Lamedtorquevirus",
    "mem": "Memtorquevirus", "samek": "Samektorquevirus",
}
BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
DNA = set("ACGTRYSWKMBDHVN")


def accession(header):
    if not header.lstrip(">").strip():
        raise ValueError("En-tête FASTA vide")
    token = header.lstrip(">").split()[0].split(",")[0]
    if "|" in token:
        parts = token.split("|")
        for i, part in enumerate(parts[:-1]):
            if part in {"gb", "ref", "emb", "dbj", "tpg", "tpe", "tpd"}:
                token = parts[i + 1]
                break
    token = re.sub(r"\.\d+$", "", token)
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", token):
        raise ValueError("Identifiant FASTA non pris en charge: " + token)
    return token


def records(text):
    result, header, seq = [], None, []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                result.append((header, "".join(seq).upper()))
            header, seq = line[1:], []
        elif header is None:
            raise ValueError("Contenu inattendu avant le premier en-tête FASTA")
        else:
            seq.append("".join(line.split()))
    if header is not None:
        result.append((header, "".join(seq).upper()))
    for header, sequence in result:
        if not sequence:
            raise ValueError("Séquence vide: " + header)
        invalid = set(sequence) - DNA
        if invalid:
            raise ValueError("Alphabet inattendu pour " + header + ": "
                             + repr("".join(sorted(invalid))))
    return result


def read_fasta(path):
    try:
        return records(Path(path).read_text())
    except ValueError as exc:
        raise ValueError(str(path) + ": " + str(exc)) from exc


def fasta_accessions(path):
    # Existing accession membership needs headers only, not sequence validation.
    # References and newly fetched sequences remain strictly validated when used.
    ids = []
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.strip()
            if line.startswith(">"):
                try:
                    ids.append(accession(line[1:]))
                except ValueError as exc:
                    raise ValueError(str(path) + ": ligne " + str(line_number)
                                     + ": " + str(exc)) from exc
            elif line and not ids:
                raise ValueError(str(path) + ": contenu avant le premier en-tête FASTA")
    if not ids:
        raise ValueError(str(path) + ": aucun en-tête FASTA")
    return ids


def write_fasta(path, rows):
    with Path(path).open("w") as handle:
        for header, sequence in rows:
            handle.write(">" + header + "\n" + sequence + "\n")


def write_tsv(path, fields, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)


def write_lines(path, lines):
    Path(path).write_text("".join(str(line) + "\n" for line in lines))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_ictv(src):
    path = src / "lookUpTable" / "ICTV.csv"
    if not path.is_file():
        raise ValueError("Table ICTV absente: " + str(path))
    statuses = {}
    with path.open(newline="") as handle:
        for number, row in enumerate(csv.reader(handle), 1):
            if not row or not any(value.strip() for value in row) or row[0].lstrip().startswith("#"):
                continue
            if len(row) != 2 or row[1].strip().lower() not in {"true", "false"}:
                raise ValueError("Ligne ICTV invalide: %s:%d" % (path, number))
            species, official = row[0].strip(), row[1].strip().lower() == "true"
            if not species or (species in statuses and statuses[species] != official):
                raise ValueError("Statut ICTV vide ou contradictoire: " + species)
            statuses[species] = official
    return statuses


def ictv_status(species, statuses):
    if "%" in species:
        return "composite"
    if species not in statuses:
        return "unknown"
    return "official" if statuses[species] else "candidate"


def classification_decision(hits, statuses):
    qualifying = {species: hit for species, hit in hits.items() if hit[0] > 69}
    if not qualifying:
        return None, "exactly_69" if any(h[0] == 69 for h in hits.values()) else "no_hit", []
    if len(qualifying) == 1:
        return next(iter(qualifying)), "single_species", []
    kinds = {species: ictv_status(species, statuses) for species in qualifying}
    if any(kind in {"unknown", "composite"} for kind in kinds.values()):
        return None, "ictv_status_review", []
    official = [s for s in qualifying if kinds[s] == "official"]
    candidates = [s for s in qualifying if kinds[s] == "candidate"]
    if not candidates:
        highest = max(hit[0] for hit in qualifying.values())
        best = [s for s, hit in qualifying.items() if hit[0] == highest]
        return (best[0], "best_ictv_identity", []) if len(best) == 1 else (None, "tied_ictv_identity", [])
    if len(official) > 1:
        return None, "candidate_multiple_ictv", []
    if official:
        return official[0], "proposed_merge_into_ictv", sorted(candidates)
    target = min(qualifying, key=lambda s: (-qualifying[s][0], s))
    return target, "proposed_candidate_merge", sorted(s for s in qualifying if s != target)


def merge_command(script_file, src, genus, sources, target):
    return shlex.join(["bash", str(script_file), "merge", genus, "--src", str(src),
                       "--from-species", *sources, "--into-species", target, "--apply"])


class NCBI:
    def __init__(self, email):
        self.email = email
        self.last_request = 0.0

    def request(self, endpoint, **params):
        params.update(tool="SCANellome_Update_2026")
        if self.email:
            params["email"] = self.email
        # Keep the key out of the command line, reports and request logs.
        if os.environ.get("NCBI_API_KEY"):
            params["api_key"] = os.environ["NCBI_API_KEY"]
        url = BASE_URL + endpoint + "?" + urllib.parse.urlencode(params)
        last_error = ""
        for attempt in range(3):
            time.sleep(max(0.0, 0.4 - (time.monotonic() - self.last_request)))
            self.last_request = time.monotonic()
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "SCANellome-2026/1.0"})
                with urllib.request.urlopen(req, timeout=45) as response:
                    return response.read().decode("utf-8")
            except urllib.error.HTTPError as exc:
                last_error = "HTTP " + str(exc.code)
                if exc.code not in {429, 500, 502, 503, 504}:
                    break
            except (urllib.error.URLError, TimeoutError, OSError) as exc:
                last_error = type(exc).__name__
            if attempt < 2:
                time.sleep(2 ** attempt)
        raise RuntimeError("NCBI inaccessible après tentative(s): " + last_error)

    def search(self, genus):
        # Same sequence-length window as the historical updater. This search
        # covers records indexed under the genus, not all unclassified primate
        # anelloviruses or candidates assigned to a different genus by SCANellome.
        term = '"' + genus + '"[Organism] AND 1800:5000[SLEN]'
        raw = self.request("esearch.fcgi", db="nuccore", term=term,
                           retmode="json", retmax=10000, idtype="acc", usehistory="y")
        payload = json.loads(raw)
        result = payload.get("esearchresult", {})
        if payload.get("error") or result.get("error") or result.get("errorlist"):
            raise RuntimeError("Recherche NCBI rejetée: " + genus)
        warning = result.get("warninglist", {})
        if warning.get("phrasesignored") or warning.get("quotedphrasesnotfound"):
            raise RuntimeError("NCBI n'a pas reconnu le terme: " + genus)
        total = int(result["count"])
        ids = list(result["idlist"])
        while len(ids) < total:
            # Use the original History set so pagination stays on one snapshot.
            page = json.loads(self.request(
                "esearch.fcgi", db="nuccore", term="#" + result["querykey"],
                WebEnv=result["webenv"], retmode="json", retmax=10000,
                retstart=len(ids), idtype="acc"))["esearchresult"]
            batch = page.get("idlist", [])
            if not batch:
                raise RuntimeError("Réponse NCBI incomplète: " + genus)
            ids.extend(batch)
        if len(ids) != total or len(set(ids)) != total:
            raise RuntimeError("Liste NCBI incomplète ou dupliquée: " + genus)
        if any(not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*[0-9](?:\.[0-9]+)?", x) for x in ids):
            raise RuntimeError("NCBI n'a pas renvoyé des accessions valides: " + genus)
        return ids, {"query": term, "query_translation": result.get("querytranslation"),
                     "count": total, "warnings": warning}

    def fetch(self, ids):
        result = []
        for start in range(0, len(ids), 100):
            batch = ids[start:start + 100]
            text = self.request("efetch.fcgi", db="nuccore", id=",".join(batch),
                                rettype="fasta", retmode="text")
            found = records(text)
            if len(found) != len(batch) or {accession(h) for h, s in found} != {accession(x) for x in batch}:
                raise RuntimeError("Téléchargement FASTA incomplet ou inattendu")
            result.extend(found)
        return result


def genus_directory(src, genus):
    candidates = [src / ("db_" + genus), src / ("db_" + genus + "-NA")]
    found = [path for path in candidates if path.is_dir()]
    if len(found) == 1:
        return found[0]
    if not found:
        raise ValueError("Dossier absent: db_" + genus + " ou db_" + genus + "-NA dans " + str(src))
    raise ValueError("Deux dossiers trouvés pour " + genus + ": "
                     + ", ".join(path.name for path in found)
                     + ". Indiquer avec --src une base contenant un seul dossier pour ce genre.")


def inventory(src):
    by_genus, elsewhere = {}, set()
    for directory in sorted(src.glob("db_*")):
        if not directory.is_dir() or directory.name.startswith("db_Group"):
            continue
        genus = directory.name[3:].split("-")[0]
        ids = by_genus.setdefault(genus, set())
        for fasta in sorted(directory.glob("*/*.fasta")):
            ids.update(fasta_accessions(fasta))
        elsewhere.update(ids)
    history = {accession(x.strip()) for x in (src / "acc.csv").read_text().splitlines() if x.strip()}
    excluded = set()
    for path in (src / "AC_list").glob("*.txt"):
        excluded.update(accession(x.strip()) for x in path.read_text().splitlines()
                        if x.strip() and not x.startswith("#"))
    return by_genus, elsewhere, history, excluded


def partition(ids, here, anywhere, history, excluded):
    groups = {key: [] for key in ("present_here", "present_elsewhere", "excluded", "history_only", "new")}
    seen = set()
    for full_id in ids:
        key = accession(full_id)
        if key in seen:
            raise ValueError("Plusieurs versions de la même accession dans l'entrée: " + key)
        seen.add(key)
        category = ("excluded" if key in excluded else "present_here" if key in here
                    else "present_elsewhere" if key in anywhere else "history_only" if key in history else "new")
        groups[category].append(full_id)
    return groups


def extract_orfs(rows, src, folder):
    work = folder / "orf_work"
    work.mkdir()
    ids = {"Q%08d" % (i + 1): accession(header) for i, (header, seq) in enumerate(rows)}
    write_fasta(work / "input.fa", [(key, rows[i][1]) for i, key in enumerate(ids)])
    (work / "ORF.fasta").touch()
    # Short temporary identifiers avoid the legacy underscore/reverse-header bug.
    with (folder / "orf.log").open("w") as log:
        for command in (["bash", str(src / "ORF1.sh"), "input.fa"],
                        ["bash", str(src / "codon_terminaison_check.sh")]):
            subprocess.run(command, cwd=work, stdout=log, stderr=subprocess.STDOUT, check=True)
    candidates = {key: set() for key in ids}
    for header, seq in records((work / "ORF_T.fasta").read_text()):
        match = re.fullmatch(r"(Q\d+)(?:_ORF1_fr[Aa]me[123])?(?:_rev)?", header)
        if match is None or match[1] not in ids:
            raise RuntimeError("En-tête ORF inattendu: " + header)
        candidates[match[1]].add(seq)
    unique, rejected, status = [], [], []
    for key, acc in ids.items():
        seqs = sorted(candidates[key])
        if len(seqs) == 1:
            unique.append((acc, seqs[0]))
        else:
            status.append((acc, "no_complete_ORF" if not seqs else "multiple_ORFs", len(seqs)))
            rejected.extend((acc + "_candidate" + str(i + 1), s) for i, s in enumerate(seqs))
    write_fasta(folder / "ORF1_candidates.fasta", unique)
    write_fasta(folder / "multiple_ORFs.fasta", rejected)
    write_tsv(folder / "orf_review.tsv", ["accession", "reason", "candidates"], status)
    return unique, status


class Comparator:
    def __init__(self, muscle, work):
        self.muscle, self.work = muscle, work
        work.mkdir()
        # Smoke test checks both the CLI and the length field used historically.
        if self.identity("ACGT" * 25, "ACGT" * 25) != 100:
            raise RuntimeError("Le test de compatibilité MUSCLE 3 a échoué")

    def identity(self, reference, query):
        inp, out = self.work / "pair.fasta", self.work / "pair.aln"
        write_fasta(inp, [("REF", reference), ("QUERY", query)])
        out.unlink(missing_ok=True)
        run = subprocess.run([self.muscle, "-in", str(inp), "-clw", "-out", str(out)],
                             capture_output=True, text=True, timeout=600)
        (self.work / "last_muscle.log").write_text(run.stdout + run.stderr)
        if run.returncode or not out.is_file():
            raise RuntimeError("MUSCLE a échoué. MUSCLE 3 (-in/-clw/-out) est requis; voir last_muscle.log")
        match = re.search(r"max length\s+(\d+)", run.stderr)
        if not match or int(match[1]) != max(len(reference), len(query)):
            raise RuntimeError("Journal MUSCLE incompatible avec la définition historique de l'identité")
        text = out.read_text()
        if not text.startswith(("MUSCLE", "CLUSTAL")):
            raise RuntimeError("Sortie CLUSTAL inattendue")
        # Identical columns / maximum raw length, without rounding.
        self.last_matches = text.count("*")
        self.last_denominator = int(match[1])
        return Fraction(self.last_matches * 100, self.last_denominator)


def classify(rows, src, genus, folder, comparator, statuses, script_file):
    base = genus_directory(src, genus)
    references = []
    for path in sorted(base.glob("*/*.fasta")):
        for header, sequence in read_fasta(path):
            references.append((path.parent.name, accession(header), sequence))
    if not references:
        raise RuntimeError("Aucune référence pour " + genus)
    accepted, pending = [], sorted(rows)
    assigned_rows, review_rows, report_rows, ambiguous_hits = [], {}, {}, {}
    cache = {}
    unknown = sorted({s for s, a, seq in references if ictv_status(s, statuses) == "unknown"})
    if unknown:
        print("  Statut ICTV inconnu (absent de lookUpTable/ICTV.csv): " + ", ".join(unknown), flush=True)
    pass_number = 0
    while pending:
        pass_number += 1
        remaining = []
        added = 0
        for acc, seq in pending:
            hits = {}
            for species, ref_acc, reference in references:
                key = (reference, seq)
                if key not in cache:
                    comparator.identity(reference, seq)
                    matches, denominator = comparator.last_matches, comparator.last_denominator
                    cache[key] = (Fraction(100 * matches, denominator), matches, denominator)
                score, matches, denominator = cache[key]
                if (species not in hits or score > hits[species][0]
                        or (score == hits[species][0] and ref_acc < hits[species][1])):
                    hits[species] = (score, ref_acc, matches, denominator)
            species, decision, sources = classification_decision(hits, statuses)
            labels = sorted(s for s, hit in hits.items() if hit[0] > 69)
            needs_review = decision not in {"single_species", "best_ictv_identity"}
            if len(labels) > 1 or decision == "exactly_69":
                relevant = labels or sorted(s for s, hit in hits.items() if hit[0] == 69)
                command = merge_command(script_file, src, genus, sources, species) if sources else ""
                report_rows[acc] = (acc, ";".join(relevant), decision, species or "",
                                    ";".join(sources), command)
                ambiguous_hits[acc] = [
                    (acc, s, hits[s][1], pass_number, hits[s][2], hits[s][3],
                     "%.8f" % hits[s][0], int(hits[s][0]), ictv_status(s, statuses)) for s in relevant]
            else:
                report_rows.pop(acc, None)
                ambiguous_hits.pop(acc, None)
            if not needs_review:
                target = folder / "additions" / base.name / species / (acc + ".fasta")
                target.parent.mkdir(parents=True, exist_ok=True)
                write_fasta(target, [(acc, seq)])
                accepted.append((acc, species, pass_number))
                score, ref_acc, matches, denominator = hits[species]
                assigned_rows.append((acc, species, pass_number, ref_acc, "%.8f" % score,
                                      matches, denominator, ictv_status(species, statuses), decision))
                references.append((species, acc, seq))
                review_rows.pop(acc, None)
                added += 1
            else:
                remaining.append((acc, seq))
                review_rows[acc] = decision
        pending = remaining
        if not added:
            break
    ambiguous = [(acc, seq, report_rows[acc][1]) for acc, seq in pending if review_rows[acc] != "no_hit"]
    pending = [(acc, seq) for acc, seq in pending if review_rows[acc] == "no_hit"]
    write_tsv(folder / "assigned.tsv", ["accession", "species", "pass", "reference", "identity_percent",
              "identical_columns", "max_raw_length", "ictv_status", "decision"], assigned_rows)
    write_fasta(folder / "no_sp.fa", pending)
    write_fasta(folder / "notOK_sp.fa", [(acc, seq) for acc, seq, labels in ambiguous])
    write_tsv(folder / "ambiguous.tsv", ["accession", "matching_species", "decision", "target_species",
              "merge_sources", "command_if_validated"], [report_rows[a] for a in sorted(report_rows)])
    write_tsv(folder / "ambiguous_hits.tsv",
              ["accession", "species", "reference", "pass", "identical_columns",
               "max_raw_length", "identity_percent", "integer_score", "ictv_status"],
              [row for a in sorted(ambiguous_hits) for row in ambiguous_hits[a]])
    return accepted, pending, ambiguous


def apply_additions(src, output, expected_acc_hash, no_complete_orfs=()):
    # Preserve history, existing FASTA accessions and completed no-ORF reviews.
    additions = sorted(output.glob("*/additions/db_*/*/*.fasta"))
    if digest(src / "acc.csv") != expected_acc_hash:
        raise RuntimeError("acc.csv a changé pendant le traitement; résultats conservés sans intégration")
    targets = [(p, src / p.parents[1].name / p.parent.name / p.name) for p in additions]
    for source, target in targets:
        if target.exists() or not target.parent.is_dir():
            raise RuntimeError("Destination existante ou dossier absent: " + str(target))
    old_acc = (src / "acc.csv").read_bytes()
    known = {accession(line.strip()) for line in old_acc.decode().splitlines() if line.strip()}
    existing = set()
    for fasta in sorted(src.glob("db_*/*/*.fasta")):
        existing.update(fasta_accessions(fasta))
    recovered = sorted(existing - known)
    accepted = {accession(records(p.read_text())[0][0]) for p in additions}
    no_complete_orfs = {accession(value) for value in no_complete_orfs}
    recorded_no_orf = sorted(no_complete_orfs - known - existing - accepted)
    missing = sorted((existing | accepted | no_complete_orfs) - known)
    if not additions and not missing:
        return 0, [], []
    (output / "acc.before.csv").write_bytes(old_acc)
    new_acc = old_acc
    if missing:
        new_acc += b"\n" if old_acc and not old_acc.endswith(b"\n") else b""
        new_acc += "".join(x + "\n" for x in missing).encode()
    created = []
    fd, temp_name = tempfile.mkstemp(prefix=".acc2026-", dir=src)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(new_acc)
            stream.flush()
            os.fsync(stream.fileno())
        for source, target in targets:
            with target.open("xb") as stream:
                created.append(target)
                stream.write(source.read_bytes())
        os.replace(temp_name, src / "acc.csv")
    except BaseException:
        for path in created:
            path.unlink(missing_ok=True)
        raise
    finally:
        Path(temp_name).unlink(missing_ok=True)
    if created:
        write_lines(output / "applied_files.txt", [p.relative_to(src) for p in created])
    return len(created), recovered, recorded_no_orf


def atomic_bytes(path, data):
    fd, temporary = tempfile.mkstemp(prefix=".update2026-", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temporary, path.stat().st_mode & 0o777)
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def merge_species(args, src, genus):
    base = genus_directory(src, genus)
    statuses = read_ictv(src)
    target = args.into_species
    sources = sorted(set(args.from_species))
    if target in sources:
        raise ValueError("La destination figure aussi parmi les sources")
    for species in sources + [target]:
        if (not re.fullmatch(re.escape(genus) + r"-[A-Za-z0-9_.-]+", species)
                or not (base / species).is_dir() or (base / species).is_symlink()):
            raise ValueError("Dossier d'espèce absent ou non pris en charge: " + species)
    if any(ictv_status(s, statuses) != "candidate" for s in sources):
        raise ValueError("Chaque source doit être explicitement False dans lookUpTable/ICTV.csv. "
                         "Une espèce officielle, composite ou de statut inconnu ne peut pas être absorbée.")
    if ictv_status(target, statuses) not in {"official", "candidate"}:
        raise ValueError("Statut ICTV de la destination inconnu")
    if base.is_symlink():
        raise ValueError("Dossier du genre symbolique non pris en charge")
    existing, filenames, additions, duplicates = {}, set(), [], []
    snapshots = {}
    for species in [target] + sources:
        paths = sorted((base / species).iterdir())
        if not paths:
            raise ValueError("Dossier d'espèce vide: " + species)
        snapshots[species] = {}
        for path in paths:
            if path.is_symlink() or not path.is_file() or path.suffix != ".fasta":
                raise ValueError("Élément inattendu dans le dossier d'espèce: " + str(path))
            rows = read_fasta(path)
            if len(rows) != 1:
                raise ValueError("La fusion demande une séquence par fichier: " + str(path))
            acc, seq = accession(rows[0][0]), rows[0][1]
            snapshots[species][path.name] = digest(path)
            if acc in existing:
                if existing[acc] != seq:
                    raise ValueError("Séquences différentes pour la même accession: " + acc)
                if species != target:
                    duplicates.append(str(path.relative_to(src)))
                continue
            if path.name in filenames:
                raise ValueError("Collision de noms de fichiers: " + path.name)
            existing[acc] = seq
            filenames.add(path.name)
            if species != target:
                additions.append((path, base / target / path.name))
    lookup = src / "lookUpTable" / "ICTV.csv"
    lookup_before = lookup.read_bytes()
    lookup_after = "".join(
        line for line in lookup_before.decode().splitlines(keepends=True)
        if not (next(csv.reader([line]), []) or [""])[0].strip() in sources).encode()
    plan = {"action": "merge", "genus": genus, "sources": sources, "target": target,
            "target_ictv": ictv_status(target, statuses), "src": str(src),
            "added_files": [str(dst.relative_to(src)) for old, dst in additions],
            "identical_accessions": duplicates, "applied": False}
    print("Fusion proposée: " + ", ".join(sources) + " -> " + target, flush=True)
    print("%d fichiers à ajouter; %d accessions déjà présentes à séquence identique." %
          (len(additions), len(duplicates)), flush=True)
    if not args.apply:
        print("Aucune modification. Ajouter --apply pour exécuter cette fusion.", flush=True)
        return 0
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    output = (args.out or src / "updates" / ("2026_merge_" + stamp + "_" + str(os.getpid()))).resolve()
    if output == base.resolve() or base.resolve() in output.parents or output == lookup.parent.resolve():
        raise ValueError("Le dossier de résultats ne doit pas être dans la base du genre ou lookUpTable")
    output.mkdir(parents=True, exist_ok=False)
    backup = output / "backup" / base.name
    backup.mkdir(parents=True)
    (output / "ICTV.before.csv").write_bytes(lookup_before)
    (output / "merge.json").write_text(json.dumps(plan, indent=2) + "\n")
    for species, expected in snapshots.items():
        paths = list((base / species).iterdir())
        if any(not p.is_file() or p.is_symlink() for p in paths):
            raise RuntimeError("Le contenu du dossier a changé: " + species)
        current = {p.name: digest(p) for p in paths}
        if current != expected:
            raise RuntimeError("La base a changé pendant la préparation: " + species)
    if lookup.read_bytes() != lookup_before:
        raise RuntimeError("ICTV.csv a changé pendant la préparation")
    created, moved, lookup_changed = [], [], False
    try:
        for source, destination in additions:
            with destination.open("xb") as handle:
                created.append(destination)
                handle.write(source.read_bytes())
        for species in sources:
            old, saved = base / species, backup / species
            os.replace(old, saved)
            moved.append((old, saved))
        atomic_bytes(lookup, lookup_after)
        lookup_changed = True
        plan["applied"] = True
        (output / "merge.json").write_text(json.dumps(plan, indent=2) + "\n")
    except BaseException:
        if lookup_changed:
            atomic_bytes(lookup, lookup_before)
        for old, saved in reversed(moved):
            os.replace(saved, old)
        for path in created:
            path.unlink(missing_ok=True)
        raise
    print("Fusion appliquée. Dossiers sources et table avant fusion conservés dans: " + str(output), flush=True)
    print("Les séquences en revue ne sont pas intégrées par merge. Relancer update --fasta sur notOK_sp.fa.", flush=True)
    return 0


def cleanup_genus(folder, applied=False):
    """Only called for completed genera; failed genera keep all diagnostics."""
    for name in ("orf_work", "align_work"):
        path = folder / name
        if path.is_dir():
            shutil.rmtree(path)
    if applied:
        path = folder / "additions"
        if path.is_dir():
            shutil.rmtree(path)
    # Headers are already present in new_sequences.fasta. Accepted ORFs are in
    # additions/ or src/; unassigned/ambiguous ORFs have their own retained FASTA.
    for name in ("orf.log", "downloaded_headers.tsv", "ORF1_candidates.fasta"):
        (folder / name).unlink(missing_ok=True)
    # Keep useful reports only; summary.tsv records zero counts explicitly.
    for name in ("multiple_ORFs.fasta", "no_sp.fa", "notOK_sp.fa",
                 "new_sequences.fasta", "orf_review.tsv", "ambiguous.tsv",
                 "ambiguous_hits.tsv", "assigned.tsv"):
        path = folder / name
        if path.is_file():
            content = path.read_text().strip()
            if not content or (path.suffix == ".tsv" and len(content.splitlines()) == 1):
                path.unlink()


def main(argv):
    script_dir = Path(argv[0]).resolve()
    script_file = Path(argv[1]).resolve()
    parser = argparse.ArgumentParser(
        prog=script_file.name,
        description="SCANellome 2026: six petits genres, comparaison à la base complète dev.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Exemples, depuis la racine Anelloviruses:
  bash src/Update-2026.sh check small
  bash src/Update-2026.sh check He Sade
  bash src/Update-2026.sh download Yod --email prenom@institut.ch
  bash src/Update-2026.sh update He --email prenom@institut.ch
  bash src/Update-2026.sh update He --fasta nouvelles.fasta --apply

Genres: He, Sade, Yod, Lamed, Mem, Samek (noms complets acceptés).
small = ces six genres. Le mode par défaut est check small.
Les dossiers db_<Genre> et db_<Genre>-NA sont acceptés; leur nom est conservé.
check: listes et résumé, sans téléchargement de séquences ni MUSCLE.
L'inventaire local lit les en-têtes FASTA; ce n'est pas un contrôle des séquences.
download: télécharge seulement les nouvelles accessions.
update: télécharge, cherche ORF1 et prépare les ajouts classés.
--apply: intègre les ajouts non ambigus et leurs accessions dans src/.
Avec update --apply, acc.csv est aussi complété avec les accessions des FASTA
déjà présents dans tous les dossiers db_* de src, tous genres confondus.
Cette synchronisation fonctionne même sans nouvelle séquence à intégrer.
Elle conserve les entrées historiques et sauvegarde acc.before.csv si modification.
update --apply inscrit aussi dans acc.csv les cas no_complete_ORF examinés avec
succès. Ils ne sont plus retentés aux lancements suivants; leurs séquences et
orf_review.tsv restent dans les résultats. Les autres cas en revue restent à traiter.
Par défaut: nettoyage des intermédiaires des genres terminés; fichiers vides retirés.
Les échecs conservent leurs diagnostics. Les ajouts restent sans --apply ou si échec.
--keep-intermediates: conserver les intermédiaires et les copies intégrées.
accessions.tsv regroupe les listes d'accessions et leurs statuts.
ambiguous_hits.tsv: meilleur hit par espèce pour les cas à plusieurs espèces
et les égalités exactes au seuil; ictv_status provient de lookUpTable/ICTV.csv.
identité précise = colonnes identiques / longueur brute maximale * 100.
Le seuil est strictement >69 %, sans troncature; integer_score est informatif.
Les égalités exactes à 69 % sans autre hit >69 % restent en revue.
Entre espèces officielles, le maximum unique détermine l'affectation.
Avec plusieurs espèces, les égalités entre maxima et les statuts inconnus
ou composites restent en revue.
ambiguous.tsv: décision, destination et commande de fusion si validation.
Les fusions candidates/candidates ou candidates/une officielle sont proposées.
Une candidate liée à plusieurs officielles reste en revue, sans commande de fusion.
Une proposition repose sur les liens de la nouvelle séquence avec ces espèces.
merge GENRE --from-species ESPÈCE... --into-species ESPÈCE: aperçu de fusion.
Ajouter --apply à merge pour exécuter la fusion; aucune fusion pendant update.
Les sources doivent être False dans ICTV.csv, la destination True ou False.
Une absence dans ICTV.csv donne unknown; les noms contenant % donnent composite.
La fusion conserve la destination, retire les lignes sources de la table et
sauvegarde les dossiers sources. acc.csv et les FASTA publiés restent inchangés.
Valider une fusion puis relancer update --fasta sur notOK_sp.fa pour reclasser.
Les pourcentages concernent les nouveaux lancements, pas les anciens résultats.

Requête: genre NCBI + longueur 1800:5000, sans filtre de date.

Extraction ORF >=1800 nt, ATG/ACG, stop terminal; alignements MUSCLE 3.
Chaque ajout classé devient une référence pour les séquences suivantes.
Les séquences non affectées sont recomparées après de nouveaux ajouts.
Les affectations déjà validées ne sont pas réattribuées pendant ce lancement.
Les ORF multiples sont mises en revue. Pas de gestion de l'origine circulaire.
Les hôtes, génotypes, nouvelles espèces, clustering et releases restent à curer.
Les fichiers FASTA publiés et les anciens scripts ne sont pas modifiés.
Python >=3.8 seul suffit pour check/download; aucune bibliothèque Python externe.
NCBI_EMAIL et NCBI_API_KEY sont lus dans l'environnement s'ils sont présents.
""")
    parser.add_argument("action", nargs="?", choices=["check", "download", "update", "merge"], default="check")
    parser.add_argument("genera", nargs="*", default=["small"])
    parser.add_argument("--src", type=Path, default=script_dir, help="dossier src du clone dev")
    parser.add_argument("--out", type=Path, help="nouveau dossier de résultats (doit être absent)")
    parser.add_argument("--email", default=os.environ.get("NCBI_EMAIL", ""))
    parser.add_argument("--fasta", type=Path, help="entrée locale, pour un seul genre; aucun appel NCBI")
    parser.add_argument("--muscle", default=os.environ.get("MUSCLE", "muscle"))
    parser.add_argument("--apply", action="store_true", help="intégrer les ajouts, ou exécuter la fusion demandée en mode merge")
    parser.add_argument("--from-species", nargs="+", help="mode merge: espèces candidates à absorber")
    parser.add_argument("--into-species", help="mode merge: espèce conservée")
    parser.add_argument("--keep-intermediates", action="store_true",
                        help="conserver les fichiers temporaires et les copies après intégration")
    args = parser.parse_args(argv[2:])
    aliases = dict(GENERA)
    aliases.update({g.lower(): g for g in GENERA.values()})
    selected = []
    for name in args.genera or ["small"]:
        values = list(GENERA.values()) if name.lower() in {"small", "all"} else [aliases.get(name.lower())]
        if None in values:
            parser.error("Genre inconnu. Choix: He Sade Yod Lamed Mem Samek ou small")
        selected.extend(g for g in values if g not in selected)
    if args.fasta and len(selected) != 1:
        parser.error("--fasta demande un seul genre")
    if args.apply and args.action not in {"update", "merge"}:
        parser.error("--apply est réservé aux modes update et merge")
    if args.action == "merge":
        if len(selected) != 1 or not args.from_species or not args.into_species or args.fasta:
            parser.error("merge demande un genre, --from-species et --into-species, sans --fasta")
    elif args.from_species or args.into_species:
        parser.error("--from-species et --into-species sont réservés au mode merge")
    src = args.src.resolve()
    if not (src / "acc.csv").is_file():
        parser.error("acc.csv absent: placer le script dans src/ ou utiliser --src /chemin/Anelloviruses/src")
    for genus in selected:
        try:
            genus_directory(src, genus)
        except ValueError as exc:
            parser.error(str(exc))
    if args.action == "merge":
        lock = os.open(src, os.O_RDONLY)
        try:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                raise RuntimeError("Une autre exécution Update-2026 utilise ce dossier src")
            return merge_species(args, src, selected[0])
        finally:
            os.close(lock)
    if args.fasta and not args.fasta.is_file():
        parser.error("FASTA local introuvable")
    if args.action == "update":
        if shutil.which(args.muscle) is None:
            parser.error("MUSCLE 3 absent; indiquer --muscle /chemin/muscle3")
        for script in ["ORF1.sh", "codon_terminaison_check.sh"]:
            if not (src / script).is_file():
                parser.error("Script manquant: " + script)
        for command in ["bash", "sed", "tr", "rev", "cut", "cat", "grep", "rm"]:
            if shutil.which(command) is None:
                parser.error("Commande manquante: " + command)
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    output = (args.out or src / "updates" / ("2026_" + args.action + "_" + stamp + "_" + str(os.getpid()))).resolve()
    output.mkdir(parents=True, exist_ok=False)
    print("Résultats: " + str(output), flush=True)
    # Directory flock avoids adding a lock file in the repository for a check.
    lock = os.open(src, os.O_RDONLY)
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        os.close(lock)
        raise RuntimeError("Une autre exécution Update-2026 utilise ce dossier src")
    try:
        before_hash = digest(src / "acc.csv")
        statuses = read_ictv(src) if args.action == "update" else {}
        ictv_hash = digest(src / "lookUpTable" / "ICTV.csv") if args.action == "update" else None
        by_genus, anywhere, history, excluded = inventory(src)
        ncbi = NCBI(args.email)
        summary, failures, no_complete_orfs = [], [], set()
        fields = ["genus", "status", "local_sequences", "ncbi_or_input_hits", "present_here",
                  "present_elsewhere", "excluded", "history_only", "new_accessions", "assigned", "review"]
        for genus in selected:
            folder = output / genus
            folder.mkdir()
            try:
                print("Recherche: " + genus, flush=True)
                if args.fasta:
                    downloaded = read_fasta(args.fasta)
                    ids = [h.split()[0] for h, s in downloaded]
                    search_info = {"source": str(args.fasta.resolve()), "sha256": digest(args.fasta)}
                else:
                    ids, search_info = ncbi.search(genus)
                    downloaded = None
                groups = partition(ids, by_genus[genus], anywhere, history, excluded)
                write_tsv(folder / "accessions.tsv", ["accession", "status"],
                          [(value, category) for category, values in groups.items() for value in values])
                if args.keep_intermediates:
                    for category, values in groups.items():
                        write_lines(folder / (category + ".accessions.txt"), values)
                (folder / "search.json").write_text(json.dumps(search_info, indent=2) + "\n")
                accepted, review_count = [], "NA"
                if args.action in {"download", "update"}:
                    needed = {accession(x) for x in groups["new"]}
                    rows = ([r for r in downloaded if accession(r[0]) in needed]
                            if downloaded is not None else ncbi.fetch(groups["new"]))
                    write_fasta(folder / "new_sequences.fasta", rows)
                    write_tsv(folder / "downloaded_headers.tsv", ["accession", "original_header"],
                              [(accession(h), h) for h, s in rows])
                    if args.action == "update":
                        if rows:
                            comparator = Comparator(args.muscle, folder / "align_work")
                            orfs, orf_review = extract_orfs(rows, src, folder)
                            accepted, pending, ambiguous = classify(orfs, src, genus, folder, comparator, statuses, script_file)
                            review_count = len(orf_review) + len(pending) + len(ambiguous)
                            no_complete_orfs.update(acc for acc, reason, _ in orf_review
                                                    if reason == "no_complete_ORF")
                        else:
                            write_tsv(folder / "assigned.tsv", ["accession", "species", "pass"], [])
                            review_count = 0
                if not args.keep_intermediates:
                    cleanup_genus(folder)
                assigned_count = len(accepted) if args.action == "update" else "NA"
                summary.append([genus, "OK", len(by_genus[genus]), len(ids),
                                *[len(groups[x]) for x in ("present_here", "present_elsewhere", "excluded", "history_only", "new")],
                                assigned_count, review_count])
                print("  %d trouvées; %d nouvelles; %d dans acc.csv seulement" %
                      (len(ids), len(groups["new"]), len(groups["history_only"])), flush=True)
            except (Exception, KeyboardInterrupt) as exc:
                message = str(exc) or type(exc).__name__
                (folder / "ERROR.txt").write_text(message + "\n")
                failures.append((genus, message))
                summary.append([genus, "ERROR", len(by_genus.get(genus, [])), *["NA"] * 8])
                print("  ERREUR: " + message, file=sys.stderr, flush=True)
                if isinstance(exc, KeyboardInterrupt):
                    break
            finally:
                write_tsv(output / "summary.tsv", fields, summary)
        manifest = {"action": args.action, "genera": selected, "src": str(src),
                    "started": stamp, "acc_sha256": before_hash, "apply_requested": args.apply,
                    "failures": failures, "applied": False,
                    "keep_intermediates": args.keep_intermediates,
                    "ictv_sha256": ictv_hash,
                    "identity": "100 * identical_columns / max_raw_length",
                    "threshold": ">69 (exact ratio)",
                    "scope": "NCBI genus-indexed records, 1800:5000 nt; accession versions ignored"}
        (output / "run.json").write_text(json.dumps(manifest, indent=2) + "\n")
        if failures:
            print("Recherche/traitement incomplet. Voir ERROR.txt. Aucun ajout intégré.", file=sys.stderr)
            return 1
        if args.apply:
            if digest(src / "lookUpTable" / "ICTV.csv") != ictv_hash:
                raise RuntimeError("ICTV.csv a changé pendant le traitement; aucun ajout intégré")
            count, recovered, recorded_no_orf = apply_additions(src, output, before_hash, no_complete_orfs)
            manifest.update(applied=True, applied_sequences=count,
                            acc_recovered_accessions=recovered,
                            acc_no_complete_orf_accessions=recorded_no_orf)
            (output / "run.json").write_text(json.dumps(manifest, indent=2) + "\n")
            if not args.keep_intermediates:
                for genus in selected:
                    cleanup_genus(output / genus, applied=True)
            print("Séquences ajoutées à src/: " + str(count), flush=True)
            print("Accessions déjà présentes ajoutées à acc.csv: " + str(len(recovered)), flush=True)
            print("Accessions sans ORF complète ajoutées à acc.csv: " + str(len(recorded_no_orf)), flush=True)
        elif args.action == "update":
            print("Ajouts préparés dans */additions/. La base reste inchangée.", flush=True)
        print("Résumé: " + str(output / "summary.tsv"), flush=True)
        return 0
    finally:
        os.close(lock)


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv[1:]))
    except (Exception, KeyboardInterrupt) as exc:
        print("ERREUR: " + (str(exc) or type(exc).__name__), file=sys.stderr)
        sys.exit(1)
PYTHON
