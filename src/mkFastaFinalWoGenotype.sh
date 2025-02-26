#! /bin/bash 
#Header ->
#GENBANK,GENUS=<>,GROUP=<>,SPECIES=<>,GENOTYPE=<>,HOST=<>,ICTV=bool

gp=$1 #group
group=$(echo $1 | cut -f2 -d "-")
ge=$(echo $1 | cut -f1 -d "-")
for species in $(ls db_$gp); do for fasta in $(ls db_$gp/$species); do for line in $(cat db_$gp/$species/$fasta) ; 
  do 
  ictv="False" #defaut
  if [ "$(grep $species"," lookUpTable/ICTV.csv | cut -d "," -f2)" = "True" ]
    then ictv="True"
    fi
  if [ ! "$(echo $fasta | cut -d "." -f2 )" = "fasta" ]
    then host=$(echo $(grep $(echo $fasta | cut -d "." -f2 )"," lookUpTable/HOST.csv) | cut -d "," -f2)
    else host="Homo_sapiens"
    fi
  if [ "${line::1}" = ">" ]
    then echo $line",GENUS="$ge",GROUP=$group,SPECIES="$species",GENOTYPE=n/a,HOST="$host",ICTV="$ictv
    else echo $line
    fi
  done; done; done
