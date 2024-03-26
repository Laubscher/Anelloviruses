import os
import sys
from Bio import Entrez

Entrez.email = "####@#####"  # Tell NCBI who you are

if sys.argv[1] == "A" :    #Alphatorquevirus
  handle = Entrez.esearch(db="nucleotide",RetMax=1000000, term='(((1800:5000[Sequence Length]) AND txid687331[Organism] ))', idtype="acc")

elif sys.argv[1] == "B" :  #Betatorquevirus
  handle = Entrez.esearch(db="nucleotide",RetMax=1000000, term='(((1800:5000[Sequence Length]) AND txid687332[Organism] ))', idtype="acc")

elif sys.argv[1] == "G" :  #Gammatorquevirus
  handle = Entrez.esearch(db="nucleotide",RetMax=1000000, term='(((1800:5000[Sequence Length]) AND txid687333[Organism] ))', idtype="acc")

elif sys.argv[1] == "H" :  #Hetorquevirus
  handle = Entrez.esearch(db="nucleotide",RetMax=1000000, term='(((1800:5000[Sequence Length]) AND txid2842729[Organism] ))', idtype="acc")

elif sys.argv[1] == "O" :  #Omegatorquevirus
  handle = Entrez.esearch(db="nucleotide",RetMax=1000000, term='(((1800:5000[Sequence Length]) AND txid2842923[Organism] ))', idtype="acc")


# searching term
#
# Anelloviridae txid687329[Organism]
#
# txid687329[Organism] not txid12618 txid687387 txid687386 txid687384 txid2012676 txid869560 txid1968861 txid2016460 txid2844800 txid2218660

# 'txid687329[Organism:exp] + host="Homo sapiens"'
#
#'((((1800:5000[Sequence Length]) AND txid687332[Organism])) NOT Laubscher[Author]) '
elif sys.argv[1] == "-1" :
  handle = Entrez.esearch(db="nucleotide",RetMax=1000000, term='(((1800:5000[Sequence Length]) AND txid687329[Organism:exp] + host="Homo sapiens") NOT txid227307[Organism] NOT txid687331[Organism] NOT txid687332[Organism] NOT txid687333[Organism] NOT txid2842729[Organism] NOT txid2842923[Organism] )', idtype="acc")


record = Entrez.read(handle)

handle.close()

file2 = open("acc2.csv", "a")  # append mode

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
    file2.write(accID.split(".")[0]+"\n")
file2.close()
