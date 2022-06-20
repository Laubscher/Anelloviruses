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
PCT=69
else
PCT=$2
fi

# name of sequences file
sequences="updating_2022"

echo "Query Sequences .."
python3 getFastaFromGB.py > $sequences

echo "Formatting Sequences .."
bash 1line.sh $sequences

echo "Seek for ORF1 .."

bash ORF1.sh "1l_"$sequences
echo "              .."
bash codon_terminaison_check.sh
rm ORF.fasta

# Process one by one avoid to recursively check %
for header in $(cat ORF_T.fasta | grep ">"); do
   grep -A1 $header ORF_T.fasta > 1by1.temp
   bash compare.sh 1by1.temp $PCT $db
   bash spAttribution output.fa $db
   rm output.fa 1by1.temp
done