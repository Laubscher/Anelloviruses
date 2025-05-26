import re

class FastaEntry:
    def __init__(self, header, sequence, genus="n/a", group="n/a", species="n/a", genotype="n/a", host="n/a", ictv="n/a"):
        self.header = header
        self.sequence = sequence
        self.genus = genus
        self.group = group
        self.species = species
        self.genotype = genotype
        self.host = host
        self.ictv = ictv
    
    def __repr__(self):
        return f"FastaEntry(header={self.header}, genus={self.genus}, group={self.group}, species={self.species}, " \
               f"genotype={self.genotype}, host={self.host}, ictv={self.ictv}, sequence_length={len(self.sequence)})"

class FastaDatabase:
    def __init__(self):
        self.entries = []
    
    def add_entry(self, entry):
        self.entries.append(entry)
    
    def query(self, key, value):
        return [entry for entry in self.entries if getattr(entry, key, None) == value]

def lire_fasta(fichier):
    db = FastaDatabase()
    with open(fichier, 'r') as f:
        header = None
        sequence = ''
        metadata = {}
        
        for line in f:
            line = line.strip()
            
            if line.startswith('>'):
                if header:
                    db.add_entry(FastaEntry(header, sequence, **metadata))
                
                # Extraction du header et des métadonnées
                parts = line[1:].split(',')
                header = parts[0]  # Premier élément = Numéro d'accès
                metadata = {key: value for key, value in 
                            (re.split(r'=', p, maxsplit=1) if '=' in p else (p, 'n/a') for p in parts[1:])}
                
                # Assurer la présence de toutes les clés
                for key in ["GENUS", "GROUP", "SPECIES", "GENOTYPE", "HOST", "ICTV"]:
                    metadata.setdefault(key, "n/a")
                
                sequence = ''
            else:
                sequence += line
        
        if header:
            db.add_entry(FastaEntry(header, sequence, **metadata))
    
    return db

