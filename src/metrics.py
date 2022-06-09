#!/usr/bin/env python3

#Purpose: make tsv with all the metrics for graphics or further analyzes, and choose the "best" virus for each virus name.
#Project: SAfia, pipeline "FeVIR"
#Usage: metrics.r prefix header black_list color
import sys

###################################
                                                         
   #       Functions:
                            
###################################
def getMedian(lst):
    if lst == []:
      return 0.0
    lst=sorted(lst)
    n=len(lst)-1 
    if n%2 == 0:
     m=n+1
     return (lst[int(n/2)]+lst[int(m/2)])/2.
    else:
     return float(lst[int(n/2)])
#####################################


##################################
#                                #
#           Argument             #
#                                #
##################################

out = sys.argv[1]  # ->prefix for output and input file name  
   
header = sys.argv[2] #file with "header" information (mainly AC, opt isolate (or strain), opt family,virus name (species for now))

black_list = sys.argv[3] #file with black listed AC -> custom threshold (because of region of homology or complexity)

color = sys.argv[4] #file with color for reports

len_read = sys.argv[5]

threshold = 3*int(len_read)-1


if threshold < 225 :
  threshold = 225

##################################
#                                #
#     get Accession Number       #
#       and reads number         #
#                                #
##################################

bga=open(out+".bga","r")

reads_dico={}
AC=set()
for line in bga:
  ACnum=line.split("@")[0]
  AC.add(ACnum)
  nbRead=line.split("@")[-1].split("\t")[0]
  reads_dico[ACnum]=nbRead
bga.close()


##################################
#                                #
#      Family  (for colors)      #
#                                #
##################################

'''Get the family from the header list:'''

family_dico={}

h=open(header,"r")
for line in h:
  fam=line.split(",")[2]
  ACnum=line.split(",")[0].split(">")[1]
  family_dico[ACnum]=fam
h.close()

##################################
#                                #
#            Colors              #
#                                #
##################################


# color by family
color_dico={}

col=open(color,"r")
for line in col:
  famColor=line.split(",")[1]
  fam=line.split(",")[0]
  color_dico[fam]=famColor
col.close()

##################################
#                                #
# Black list                     #
# Virus with homology region     #
# or contamination               #
#                                #
##################################

#
regions_dico={}

b=open(black_list,"r")
for line in b:
  length=int(line.split(",")[1])
  ACnum=line.split(",")[0]
  regions_dico[ACnum]=length
b.close()

##################################
#                                #
#    Virus name & species        #
#                                #
##################################

virus_name={}
virus_sp={}
h=open(header,"r")
for line in h:
  line=line.strip()
  name=line.split(",")[1]      #virus names are species by default, but can be (clinical) usal name or group of species (eg.anellovirus)
  acess=line.split(",")[0].split(">")[1]
  if name in virus_name:
    virus_name[name].append(acess)
  else:
    virus_name[name]=[acess]
  virus_sp[acess]= line.split(",")[-1] #virus species
h.close()

##################################
#                                #
#      Coverage and Co.          #
#                                #
##################################
'''
   For each segment detected, compute the genome coverage, the depth (+ segment length, reads, color )
'''

value_Dico={}

for segment in AC :    #segment = accesion numbr. of the segment
 try: 
  cov_Dico={}

  bga=open("{}.bga".format(out),"r")
  
  for line in bga:
      if line.split("@")[0] == segment:                        # if the same name will make a dictionary of the coverage {position : coverage}
          for i in range(int(line.split("\t")[1]), int(line.split("\t")[2]), 1):
            cov_Dico[i+1]=(int(line.split("\t")[3]))

  bga.close()
  
  notZero=[v for v in cov_Dico.values() if v > 0]     #get a list with only th coverage vale > 0
  coverage= len(notZero)                
  depth = getMedian(notZero)                          #median of the depth of the covered position. 
  genomeLen=-1
  #if bool(cov_Dico):
  genomeLen = max(cov_Dico.keys())
  if family_dico[segment] in color_dico:
   colors= color_dico[family_dico[segment]]
  else:
   colors="gray50"           # default color 
  
  value_Dico[segment] = [coverage,depth,reads_dico[segment],round((coverage/float(genomeLen))*100,3),genomeLen,segment,colors]  #family_dico[segment] can be add family
 except:
  pass
##################################
#                                #
#   All and Best hit by virus    #
#                                #
##################################
allHits=[]
bestHits=[]
for key in virus_name.keys():                         #  itaration by segment name 
  allHits_forthekey=[]
  for segment in virus_name[key]:                     # (segment = AC, iterate all the Ac for one virus name)
   try:
    if segment in AC:
      value_Dico[segment].insert(0,key)               # insert virus name in fisrt pos of the value_dico
      value_Dico[segment].insert(-1,virus_sp[segment])    # insert virus species in last pos of the value dico

      if segment.split(":")[0] in regions_dico.keys():                     # test if the ACnum is in the list "FP" potential, if yes check the coverage threshold
          if value_Dico[segment][1] > regions_dico[segment.split(":")[0]]:
            allHits.append(value_Dico[segment])
            allHits_forthekey.append(value_Dico[segment])
      else:
          allHits.append(value_Dico[segment])
          allHits_forthekey.append(value_Dico[segment])

   except:
     pass
  if allHits_forthekey != [] :



    maxCov=max([item[1] for item in allHits_forthekey])   # 1 = coverage /!\ 
    allMaxCov=[]
    allMaxDepth=[]
    allMaxReads=[]
    for allV in allHits_forthekey :                       # keep the one(s) with the maximum coverage.
      if allV[1] == maxCov:# 1 = coverage
        if allV[1]>threshold:                             # minimum coverage
          allMaxCov.append(allV)
    if len(allMaxCov)==1:                                 # if only one with this coverage                                                 
        bestHits.append(allMaxCov.pop())
    elif len(allMaxCov)>1:
      maxDepth=max([item[2] for item in allMaxCov])       # keep the one(s) with the maximum depth.
      for allV in allMaxCov :                       
        if allV[2] == maxDepth:
            allMaxDepth.append(allV)
      if len(allMaxDepth) == 1:                             # if only one with this Depth
          bestHits.append(allMaxDepth.pop())
      elif len(allMaxDepth) > 1:
        maxReads = max([item[3] for item in allMaxDepth])   # keep the one(s) with the maximum reads
        for allV in allMaxDepth:
          if allV[3] == maxReads:
              allMaxReads.append(allV)
        bestHits.append(allMaxReads[0])                                      # keep the first with max Reads (avoid tie)
bestOut=open(out+"_best.tsv","w")

bestOut.write("Virus\t% Segment coverage\tDepth(Median)\tReads\tCoverage\tSegment length\tSpecies\tAcession number\tColors\n")

for i in bestHits:

 bestOut.write(str(i[0])+"\t"+str(i[4])+"\t"+str(i[2])+"\t"+str(i[3])+"\t"+str(i[1])+"\t"+str(i[5])+"\t"+str(i[7])+"\t"+str(i[6].split(":")[0])+"\t"+str(i[-1])+"\n")

bestOut.close()
allOut=open(out+"_all.tsv","w")

allOut.write("Virus\t% Segment coverage\tDepth(Median)\tReads\tCoverage\tSegment length\tSpecies\tAcession number\n")

for i in allHits:

 allOut.write(str(i[0])+"\t"+str(i[4])+"\t"+str(i[2])+"\t"+str(i[3])+"\t"+str(i[1])+"\t"+str(i[5])+"\t"+str(i[7])+"\t"+str(i[6])+"\n")


allOut.close()
