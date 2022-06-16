#! /bin/bash
#check that ORF seq finish by STOP codon
grep -s -B1 tag$ ORF.fasta | sed /--/d | tr "atgc" "ATGC" > ORF_T.fasta

grep -s -B1 tga$ ORF.fasta | sed /--/d | tr "atgc" "ATGC" >> ORF_T.fasta

grep -s -B1 taa$ ORF.fasta | sed /--/d | tr "atgc" "ATGC" >> ORF_T.fasta

#correction header name
sed -i "s/A_/A/g" ORF_T.fasta
sed -i "s/T_/T/g" ORF_T.fasta
sed -i "s/G_/G/g" ORF_T.fasta
sed -i "s/C_/C/g" ORF_T.fasta
sed -i "s/N_/N/g" ORF_T.fasta
sed -i "s/^NC/NC_/g" ORF_T.fasta