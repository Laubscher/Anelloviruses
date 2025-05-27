
def lire_fasta(fichier):
    sequences = {}
    with open(fichier, 'r') as f:
        header = None
        sequence = ''
        
        for line in f:
            line = line.strip()
            
            # Si la ligne commence par '>', c'est un en-tête de séquence
            if line.startswith('>'):
                # Si une séquence a déjà été lue, on l'ajoute au dictionnaire
                if header:
                    sequences[header] = sequence
                
                header = line[1:]  # Retirer le '>' de l'en-tête
                sequence = ''
            else:
                sequence += line
        
        # Ajouter la dernière séquence au dictionnaire
        if header:
            sequences[header] = sequence
    
    return sequences




from Bio import Align


ref=lire_fasta("ref.fasta")
all=lire_fasta("all.fasta")

for seq in all.keys():
 seqScore=0
 query=all[seq]
 for seq2 in ref.keys():
    target=ref[seq2]
  

    alignment=Align.PairwiseAligner(scoring="megablast").align(target, query)[0]

    pct=alignment.counts()[1]/len(alignment[0])
    #print(pct)
    if pct >= 0.91:
      seqScore=1 #deja rep
 if seqScore==0 :
   print(">"+seq)
   print(all[seq])


