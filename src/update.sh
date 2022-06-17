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

