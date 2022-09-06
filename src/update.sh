#!/bin/bash

if [ -z "$*" ]; then
  echo "Arguments are missing"
  exit 0
fi

if [ -z "$1" ]; then
db="db_Alphatorquevirus"                          #db_Betatorquevirus
else
db=$1
fi

if [ -z "$2" ]; then
PCT=68
else
PCT=$2
fi

if [ -z "$3" ]; then
ORF=FALSE
else
ORF=TRUE
fi

# name of sequences file
sequences="updating_2022"

#echo "Query Sequences .."
#python3 getFastaFromGB.py > $sequences

if [ "$ORF" == "FALSE" ]; then

echo "Formatting Sequences .."
bash 1line.sh $sequences

echo "Seek for ORF1 .."

bash ORF1.sh "1l_"$sequences
echo "              .."
bash codon_terminaison_check.sh
rm ORF.fasta

else cp updating_2022 ORF_T.fasta
fi
rm -f no_sp.fa
# Process one by one avoid to many recursive check %
for header in $(cat ORF_T.fasta | grep ">"); do
   grep -A1 $header ORF_T.fasta > 1by1.temp
   bash compare.sh 1by1.temp $PCT $db
   bash spAttribution.sh output.fa $db
   rm output.fa 1by1.temp
done

#recursive check until no_sp.fa stop changing

check2=$(wc -l ORF_T.fasta | cut -f1 -d" " )
check1=$(wc -l no_sp.fa | cut -f1 -d" " )

rm ORF_T.fasta

notSame=1

if [ $check1 -eq $check2 ]
  then
  notSame=0
fi

while [ $notSame -eq 1 ]
do
  check2=$(echo $check1)
   mv no_sp.fa ORF_T.fasta

   for header in $(cat ORF_T.fasta | grep ">"); do
     grep -A1 $header ORF_T.fasta > 1by1.temp
     bash compare.sh 1by1.temp $PCT $db
     bash spAttribution.sh output.fa $db
     rm output.fa 1by1.temp
   done

  check1=$(wc -l no_sp.fa | cut -f1 -d" " )
  if [ $check1 -eq $check2 ]
    then
     notSame=0
  fi
done