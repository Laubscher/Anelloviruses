#!/usr/bin/env python3

import sys
genotypeDico=dict()
fasta=open(str(sys.argv[1]), "r")
genotype=""
AC=[]
for l in fasta:
  if l[0]==">":
    genotype=set(genotype)
    try:
      genotype.remove('n/a')
    except:
      pass
    if len(genotype)==1:
      genotype=genotype.pop()
    else:
      genotype=str(genotype)
    if genotype == "set()":
      genotype='n/a'
    for i in AC:
      genotypeDico[i]=genotype
    genotype=[]
    AC=[]
  else: 
    AC.append(l.split(",")[1].split(">")[1])
    genotype.append(l.split("GENOTYPE=")[1].split(",")[0])
    
genotype=set(genotype)
try:
      genotype.remove('n/a')
except:
      pass
if len(genotype)==1:
      genotype=genotype.pop()
else:
      genotype=str(genotype)
if genotype == "set()":
      genotype='n/a'
for i in AC:
      genotypeDico[i]=genotype


fasta.close()



fasta2=open(str(sys.argv[2]), "r")
for li in fasta2:
  if li[0]==">":
    ACC=li.split(">")[1].split(",")[0]
    li2=li.split("GENOTYPE=")[0]+"GENOTYPE="+str(genotypeDico[ACC])+","+li.rstrip().split("GENOTYPE=")[1].split(",")[1]
    print(li2)
  else:
    print(li.rstrip())
fasta2.close()

