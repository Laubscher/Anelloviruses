#! /bin/bash 
#Header ->
#GENBANK,GENUS=<>,GROUP=<>,SPECIES=<>,GENOTYPE=<>,HOST=<>,ICTV=bool

gp=$1 #group
group=$(echo $1 | cut -f2 -d "-")
ge=$(echo $1 | cut -f1 -d "-")
for species in $(ls db_$gp); do for fasta in $(ls db_$gp/$species); do for line in $(cat db_$gp/$species/$fasta) ; 
  do 
  ICTV="False" #defaut
  if [ "$(grep $species"," lookUpTable/ICTV.csv | cut -d "," -f2)" = "True" ]
    then ICTV="True"
    fi
  if [ "${line::1}" = ">" ]
  then if [ "$(echo $fasta | cut -d "-" -f2 )" = "Pt.fasta" ]   
    then echo $line",GENUS="$ge",GROUP=$group,SPECIES="$species",GENOTYPE=n/a,HOST=Pan_troglodytes,ICTV="$ICTV
    else echo $line",GENUS="$ge",GROUP=$group,SPECIES="$species",GENOTYPE=n/a,HOST=Homo_sapiens,ICTV="$ICTV    #as defaut Homo_sapiens | ( n/a )
    fi
  else echo $line
  fi
  done; done; done
