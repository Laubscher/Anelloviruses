import os
#from Bio import SeqIO
from Bio import Entrez

Entrez.email = "####@#####"  # Tell NCBI who you are

handle = Entrez.esearch(db="nucleotide",RetMax=1000000, term='txid687329[Organism:exp] + host="Homo sapiens" ', idtype="acc")

record = Entrez.read(handle)

handle.close()

with open('acc.csv') as f:
  DBaccIDlist = f.read().splitlines()

for accID in record["IdList"]:
  
  if accID.split(".")[0] not in DBaccIDlist:
    #handle = Entrez.efetch(db="nucleotide", id="AB018887.1", rettype="gb", retmode="text")
    #gb=handle.read()
    #taxo=gb.split("ORGANISM")[1].split("REFERENCE")[0].split(";")[-1].rstrip() #no need for taxa
    
    handle = Entrez.efetch(db="nucleotide", id=accID, rettype="fasta", retmode="text")
    fasta=handle.read()
    print(fasta)