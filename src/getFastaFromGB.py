import os
#from Bio import SeqIO
from Bio import Entrez

Entrez.email = "####@#####"  # Tell NCBI who you are

# searching term
#
# Anelloviridae txid687329[Organism]
#
# txid687329[Organism] not txid12618 txid687387 txid687386 txid687384 txid2012676 txid869560 txid1968861 txid2016460 txid2844800 txid2218660
# -
#
# Alphatorquevirus txid687331[Organism]
#
# Betatorquevirus txid687332[Organism]
#
# Gammatorquevirus txid687333[Organism]
#
# Hetorquevirus txid????????[Organism]
#
# Omegatorquevirus txid2842923[Organism]
#
# 'txid687329[Organism:exp] + host="Homo sapiens"'
#
# -txid12618 txid687387 txid687386 txid687384 txid2012676 txid869560 txid1968861 txid2016460
# txid2844800 txid2218660

handle = Entrez.esearch(db="nucleotide",RetMax=1000000, term='((((1800:5000[Sequence Length]) AND txid687332[Organism])) NOT Laubscher[Author]) ', idtype="acc")

record = Entrez.read(handle)

handle.close()

with open('acc.csv') as f:
  DBaccIDlist = f.read().splitlines()

for accID in record["IdList"]:
  
  if accID.split(".")[0] not in DBaccIDlist:
    '''if taxa needed'''
    #handle = Entrez.efetch(db="nucleotide", id=accID, rettype="gb", retmode="text")
    #gb=handle.read()
    #taxo=gb.split("ORGANISM")[1].split("REFERENCE")[0].split(";")[-1].rstrip()
    
    handle = Entrez.efetch(db="nucleotide", id=accID, rettype="fasta", retmode="text")
    fasta=handle.read()
    print(fasta)