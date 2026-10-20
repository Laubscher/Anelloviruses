#!/usr/bin/env python3
"""Copier les nouveaux FASTA de dev et créer une prérelease par genre."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request


BRANCH = "db_Anellovirus"
PAGES_BRANCH = "gh-pages"
NAME = re.compile(r"FASTA/([A-Z][A-Za-z]+)_(20\d{2}\.\d+)\.fasta\Z")
FIELDS = {"GENUS", "GROUP", "SPECIES", "GENOTYPE", "HOST", "ICTV"}


def run(*args, cwd=None):
    result = subprocess.run(args, cwd=cwd, capture_output=True, check=False)
    if result.returncode:
        raise RuntimeError(result.stderr.decode(errors="replace").strip())
    return result.stdout


def git(*args, cwd=None):
    return run("git", *args, cwd=cwd)


def version(value):
    return tuple(map(int, value.split(".")))


def blob(ref, path):
    entry = git("ls-tree", ref, "--", path).decode().strip()
    if not entry:
        return None
    if not entry.startswith(("100644 blob ", "100755 blob ")):
        raise ValueError(f"Fichier ordinaire attendu : {path}")
    return git("show", f"{ref}:{path}")


def read_fasta(data, genus, strict=True):
    records = {}
    accession = None
    sequence = []
    metadata = {}

    def finish():
        if accession is None:
            return
        seq = "".join(sequence).upper()
        if not seq or not re.fullmatch(r"[ACGTRYSWKMBDHVN]+", seq):
            raise ValueError(f"Séquence vide ou alphabet inattendu : {accession}")
        records[accession] = (seq, metadata.copy())

    for line in data.decode("utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish()
            parts = line[1:].split(",")
            accession = parts[0]
            if not re.fullmatch(r"[A-Za-z0-9_]+(?:\.\d+)?", accession):
                raise ValueError(f"Accession invalide : {accession}")
            if accession in records:
                raise ValueError(f"Accession répétée : {accession}")
            metadata = {}
            for field in parts[1:]:
                key, separator, value = field.partition("=")
                if not separator or not value or key in metadata:
                    raise ValueError(f"En-tête invalide : {accession}")
                metadata[key] = value
            if metadata.get("GENUS") != genus:
                raise ValueError(f"GENUS différent du nom du fichier : {accession}")
            if strict and (set(metadata) != FIELDS or metadata["ICTV"] not in {"True", "False"}):
                raise ValueError(f"Champs SCANellome incomplets ou invalides : {accession}")
            sequence = []
        else:
            if accession is None:
                raise ValueError("Séquence présente avant le premier en-tête FASTA")
            sequence.append(line)
    finish()
    if not records:
        raise ValueError("FASTA vide")
    return records


def candidates(before, after):
    for sha in (before, after):
        if not re.fullmatch(r"[0-9a-f]{40}", sha) or sha == "0" * 40:
            raise ValueError("Deux commits existants sont nécessaires ; création de branche non prise en charge.")
        git("cat-file", "-e", f"{sha}^{{commit}}")
    paths = git("diff", "--name-only", "--diff-filter=A", "--no-renames", "-z", before, after, "--", "FASTA/")
    found = []
    for raw in paths.split(b"\0"):
        if not raw:
            continue
        path = raw.decode("utf-8")
        match = NAME.fullmatch(path)
        if not match or match[1] == "Anellovirus":
            continue
        genus, number = match.groups()
        data = blob(after, path)
        records = read_fasta(data, genus)
        found.append(dict(path=path, genus=genus, version=number, tag=Path(path).stem,
                          data=data, records=records))
    return sorted(found, key=lambda item: (item["genus"], version(item["version"])))


def release_notes(item, target):
    older = []
    for path in git("ls-tree", "-r", "--name-only", target, "--", "FASTA/").decode().splitlines():
        match = NAME.fullmatch(path)
        if match and match[1] == item["genus"] and version(match[2]) < version(item["version"]):
            older.append((version(match[2]), path))
    current = item["records"]
    lines = ["Changelog:", "", f'[{item["version"]}]:', ""]
    if not older:
        lines.append(f"- {len(current)} representative sequences (first published version for this genus).")
    else:
        _, previous_path = max(older)
        old = read_fasta(blob(target, previous_path), item["genus"], strict=False)
        added, removed = current.keys() - old.keys(), old.keys() - current.keys()
        shared = current.keys() & old.keys()
        lines.append(f"- +{len(added)} representative sequences ({len(current)} total).")
        if removed:
            lines.append(f"- Removed {len(removed)} previous representatives.")
        changed = sum(current[acc][0] != old[acc][0] for acc in shared)
        reassigned = sum(current[acc][1].get("SPECIES") != old[acc][1].get("SPECIES") for acc in shared)
        annotations = sum(current[acc][1] != old[acc][1] for acc in shared)
        if changed:
            lines.append(f"- Updated nucleotide sequences for {changed} existing representatives.")
        if reassigned:
            lines.append(f"- Updated species assignments for {reassigned} existing representatives.")
        ictv_only = annotations and all(
            {**old[acc][1], "ICTV": "True"} == current[acc][1] for acc in shared
        )
        if ictv_only:
            lines.append("- Added `ICTV=True` annotations to FASTA headers.")
        elif annotations:
            lines.append(f"- Updated FASTA annotations for {annotations} existing representatives.")
        lines.extend(["", f"Previous FASTA: `{Path(previous_path).name}`."])
    return "\n".join(lines) + "\n"


class GitHub:
    def __init__(self, repository, token):
        self.repository = repository
        self.token = token
        self.base = f"https://api.github.com/repos/{repository}"

    def request(self, path, method="GET", payload=None, missing=False):
        headers = {"Authorization": f"Bearer {self.token}", "Accept": "application/vnd.github+json",
                   "X-GitHub-Api-Version": "2026-03-10", "User-Agent": "SCANellome-prereleases"}
        data = None
        if payload is not None:
            data = json.dumps(payload).encode()
            headers["Content-Type"] = "application/json"
        request = urllib.request.Request(self.base + path, data=data, headers=headers, method=method)
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
                return json.load(response)
        except urllib.error.HTTPError as error:
            if missing and error.code == 404:
                return None
            raise RuntimeError(f"GitHub API : HTTP {error.code} pour {path}") from error

    def get_release(self, tag):
        release = self.request("/releases/tags/" + urllib.parse.quote(tag, safe=""), missing=True)
        if release:
            return release
        # Les brouillons sont recherchés aussi pour reprendre un envoi interrompu.
        page = 1
        while True:
            releases = self.request(f"/releases?per_page=100&page={page}")
            for release in releases:
                if release["tag_name"] == tag:
                    return release
            if len(releases) < 100:
                return None
            page += 1

    def create_release(self, item, notes):
        return self.request("/releases", "POST", {
            "tag_name": item["tag"], "target_commitish": BRANCH, "name": item["tag"], "body": notes,
            "draft": True, "prerelease": True, "make_latest": "false",
        })

    def upload(self, release, item):
        # gh gère l'authentification et l'envoi de la pièce jointe sans écrasement.
        with tempfile.TemporaryDirectory(prefix="scanellome-asset-") as temporary:
            path = Path(temporary) / Path(item["path"]).name
            path.write_bytes(item["data"])
            run("gh", "release", "upload", item["tag"], str(path), "--repo", self.repository)

    def verify_asset(self, release, item):
        name = Path(item["path"]).name
        assets = [asset for asset in release["assets"] if asset["name"] == name]
        if not assets:
            return False
        if len(assets) != 1 or assets[0]["state"] != "uploaded":
            raise ValueError(f"Pièce jointe incomplète ou répétée : {name}")
        asset = assets[0]
        expected = "sha256:" + hashlib.sha256(item["data"]).hexdigest()
        if asset.get("digest"):
            matches = asset["digest"] == expected
        else:
            with tempfile.TemporaryDirectory(prefix="scanellome-verify-") as temporary:
                run("gh", "release", "download", item["tag"], "--repo", self.repository,
                    "--pattern", name, "--dir", temporary)
                matches = (Path(temporary) / name).read_bytes() == item["data"]
        if not matches or asset["size"] != len(item["data"]):
            raise ValueError(f"Pièce jointe différente du FASTA source : {name}")
        return True

    def publish(self, release):
        return self.request(f'/releases/{release["id"]}', "PATCH", {
            "draft": False, "prerelease": True, "make_latest": "false",
        })


def tag_commit(tag):
    result = subprocess.run(["git", "show-ref", "--verify", "--quiet", f"refs/tags/{tag}"])
    if result.returncode == 1:
        return None
    if result.returncode:
        raise RuntimeError(f"Impossible de lire le tag {tag}")
    return git("rev-parse", f"refs/tags/{tag}^{{commit}}").decode().strip()


def check_version(item, target):
    versions = []
    for path in git("ls-tree", "-r", "--name-only", target, "--", "FASTA/").decode().splitlines():
        match = NAME.fullmatch(path)
        if match and match[1] == item["genus"] and path != item["path"]:
            versions.append((version(match[2]), match[2]))
    for tag in git("tag", "--list", item["genus"] + "_20*").decode().splitlines():
        match = NAME.fullmatch("FASTA/" + tag + ".fasta")
        if match:
            versions.append((version(match[2]), match[2]))
    if versions:
        latest, label = max(versions)
        if version(item["version"]) <= latest:
            raise ValueError(f'{item["tag"]} : la version doit être strictement supérieure à {label}.')


def check_existing(item, target):
    present = blob(target, item["path"])
    if present is not None and present != item["data"]:
        raise ValueError(f'{item["path"]} existe déjà avec un contenu différent. Utiliser une nouvelle version.')
    commit = tag_commit(item["tag"])
    if commit:
        if blob(commit, item["path"]) != item["data"]:
            raise ValueError(f'Tag existant incompatible : {item["tag"]}')
        git("merge-base", "--is-ancestor", commit, target)
    else:
        check_version(item, target)
    return commit


def publish_items(items, after, api, dry_run=False, published=None):
    git("fetch", "origin", f"+refs/heads/{BRANCH}:refs/remotes/origin/{BRANCH}", "--tags")
    target = f"refs/remotes/origin/{BRANCH}"
    # Vérifier tous les fichiers avant la première publication.
    for item in items:
        check_existing(item, target)
        release_notes(item, target)
    if dry_run:
        for item in items:
            print(f'À publier : {item["tag"]}\n{release_notes(item, target)}')
        return
    with tempfile.TemporaryDirectory(prefix="scanellome-publication-") as temporary:
        checkout = Path(temporary) / "publication"
        git("worktree", "add", "--detach", str(checkout), target)
        try:
            for item in items:
                current = git("rev-parse", "HEAD", cwd=checkout).decode().strip()
                commit = check_existing(item, current)
                notes = release_notes(item, f"{commit}^" if commit else current)
                release = api.get_release(item["tag"])
                if release and not commit:
                    raise ValueError(f'Prérelease existante sans tag compatible : {item["tag"]}')
                if not commit:
                    path = checkout / item["path"]
                    if blob(current, item["path"]) is None:
                        path.parent.mkdir(parents=True, exist_ok=True)
                        path.write_bytes(item["data"])
                        git("add", "--", item["path"], cwd=checkout)
                    for old_path in git("ls-tree", "-r", "--name-only", current, "--", "FASTA/").decode().splitlines():
                        match = NAME.fullmatch(old_path)
                        if match and match[1] == item["genus"] and old_path != item["path"]:
                            git("rm", "--", old_path, cwd=checkout)
                    if git("diff", "--cached", "--name-only", cwd=checkout):
                        message = (f'Add {item["genus"]} {item["version"]}: '
                                   f'{len(item["records"])} representatives\n\nSource-dev: {after}')
                        git("-c", "user.name=github-actions[bot]", "-c",
                            "user.email=41898282+github-actions[bot]@users.noreply.github.com",
                            "commit", "-m", message, cwd=checkout)
                    commit = git("rev-parse", "HEAD", cwd=checkout).decode().strip()
                    # Le commit et le tag sont poussés ensemble ; aucun push forcé.
                    git("push", "--atomic", "origin", f"{commit}:refs/heads/{BRANCH}",
                        f'{commit}:refs/tags/{item["tag"]}')
                if release is None:
                    release = api.create_release(item, notes)
                if not api.verify_asset(release, item):
                    if not release["draft"]:
                        raise ValueError(f'Release déjà publiée sans le FASTA attendu : {item["tag"]}')
                    api.upload(release, item)
                    release = api.request(f'/releases/{release["id"]}')
                    if not api.verify_asset(release, item):
                        raise ValueError(f'Pièce jointe absente après envoi : {item["tag"]}')
                if release["draft"]:
                    release = api.publish(release)
                if published is not None:
                    published.append(item)
                status = "prérelease" if release["prerelease"] else "release déjà validée"
                print(f'{item["tag"]} : {status}, {len(item["records"])} représentants. {release["html_url"]}')
                summary = os.environ.get("GITHUB_STEP_SUMMARY")
                if summary:
                    with open(summary, "a", encoding="utf-8") as handle:
                        handle.write(f'- [{item["tag"]}]({release["html_url"]}) : '
                                     f'{len(item["records"])} représentants, {status}.\n')
        finally:
            git("worktree", "remove", "--force", str(checkout))


def update_index(content, items):
    for item in items:
        genus = re.escape(item["genus"])
        pattern = re.compile(
            r'(?P<start><a\b[^>]*?\bhref=["\'])'
            r'https://github\.com/Laubscher/Anelloviruses/releases/tag/'
            rf'(?P<old>{genus}(?:_20\d{{2}}\.\d+)?)(?P<end>["\'][^>]*>)'
            r'(?P<label>[^<]*)(?P<close></a>)', re.IGNORECASE)
        matches = list(pattern.finditer(content))
        if not matches:
            print(f'::warning::Aucun lien existant pour {item["genus"]} dans gh-pages/index.html.')
            continue
        if len(matches) != 1:
            raise ValueError(f'Plusieurs liens de téléchargement pour {item["genus"]} dans index.html')
        match = matches[0]
        old = match["old"]
        old_version = old.removeprefix(item["genus"] + "_")
        # Une reprise ancienne ne doit jamais faire reculer le lien du site.
        if old_version != old and version(old_version) > version(item["version"]):
            continue
        if match["label"].count(old) != 1:
            raise ValueError(f'Texte du lien inattendu pour {item["genus"]} dans index.html')
        replacement = (match["start"] + 'https://github.com/Laubscher/Anelloviruses/releases/tag/'
                       + item["tag"] + match["end"] + match["label"].replace(old, item["tag"])
                       + match["close"])
        content = content[:match.start()] + replacement + content[match.end():]
    return content


def update_pages(items, api):
    settings = api.request("/pages")
    source = settings.get("source", {})
    if settings.get("build_type", "legacy") != "legacy" or source.get("branch") != PAGES_BRANCH or source.get("path") != "/":
        raise ValueError("GitHub Pages doit publier depuis gh-pages, dossier / (root). Les préreleases sont conservées.")
    git("fetch", "origin", f"+refs/heads/{PAGES_BRANCH}:refs/remotes/origin/{PAGES_BRANCH}")
    with tempfile.TemporaryDirectory(prefix="scanellome-pages-") as temporary:
        checkout = Path(temporary) / "pages"
        git("worktree", "add", "--detach", str(checkout), f"origin/{PAGES_BRANCH}")
        try:
            path = checkout / "index.html"
            original = path.read_text(encoding="utf-8")
            updated = update_index(original, items)
            if original != updated:
                path.write_text(updated, encoding="utf-8")
                git("add", "--", "index.html", cwd=checkout)
                git("-c", "user.name=github-actions[bot]", "-c",
                    "user.email=41898282+github-actions[bot]@users.noreply.github.com", "commit", "-m",
                    "Update database download links: " + ", ".join(item["tag"] for item in items), cwd=checkout)
                commit = git("rev-parse", "HEAD", cwd=checkout).decode().strip()
                git("push", "origin", f"{commit}:refs/heads/{PAGES_BRANCH}")
            expected = git("rev-parse", "HEAD", cwd=checkout).decode().strip()
        finally:
            git("worktree", "remove", "--force", str(checkout))
    build = api.request("/pages/builds/latest", missing=True)
    if build and build.get("commit") == expected and build.get("status") == "built":
        print("GitHub Pages est déjà à jour.")
        return
    # Un commit poussé avec GITHUB_TOKEN ne déclenche pas Pages à lui seul.
    api.request("/pages/builds", "POST")
    for _ in range(60):
        build = api.request("/pages/builds/latest", missing=True)
        if build and build.get("commit") == expected:
            if build.get("status") == "built":
                print("Index publié sur GitHub Pages : " + settings.get("html_url", ""))
                return
            if build.get("status") == "errored":
                raise RuntimeError("Échec du déploiement GitHub Pages : " + str(build.get("error")))
        time.sleep(10)
    raise RuntimeError("Déploiement GitHub Pages non confirmé après 10 minutes. Voir Actions / Pages.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", required=True)
    parser.add_argument("--after", required=True)
    parser.add_argument("--dry-run", action="store_true", help="Vérifier et afficher les notes sans publier")
    args = parser.parse_args()
    items = candidates(args.before, args.after)
    if not items:
        print("Aucun nouveau FASTA de genre à publier dans FASTA/.")
        return
    api = None
    if not args.dry_run:
        if os.environ.get("GITHUB_REF") != "refs/heads/dev" or os.environ.get("GITHUB_REPOSITORY") != "Laubscher/Anelloviruses":
            raise ValueError("Publication autorisée uniquement depuis le workflow de dev dans Laubscher/Anelloviruses.")
        api = GitHub(os.environ["GITHUB_REPOSITORY"], os.environ["GH_TOKEN"])
    published = []
    try:
        publish_items(items, args.after, api, args.dry_run, published)
    finally:
        if published:
            update_pages(published, api)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, RuntimeError, OSError) as error:
        print(f"ERREUR : {error}", file=sys.stderr)
        sys.exit(1)
