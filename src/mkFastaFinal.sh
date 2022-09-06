#! /bin/bash 
#Header ->
#GENBANK,GENUS=<>,GROUP=<>,SPECIES=<>,GENOTYPE=<>,HOST=<>

sp=$1 #species as in db_species eg. "TTV1"
ge=$2 #Genus

for genotype in $(ls db_$sp); do for fasta in $(ls db_$sp/$genotype); do for line in $(cat db_$sp/$genotype/$fasta) ; 
do 
if [ "${line::1}" = ">" ]
then echo $line",GENUS="$ge",GROUP=NA,SPECIES="$sp",GENOTYPE="$genotype",HOST=Homo_sapiens"
else echo $line
fi
done; done; done