#!/bin/bash
# name of sequences file
sequences="human_alphatorquevirus_2022_2"
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
   bash compare.sh 1by1.temp
   bash spAttribution output.fa
   rm output.fa 1by1.temp
done
