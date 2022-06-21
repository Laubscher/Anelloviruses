#! /bin/bash
#check that ORF seq finish by STOP codon
grep -s -B1 tag$ ORF.fasta | sed /--/d | tr "atgc" "ATGC" > ORF_T.fasta

grep -s -B1 tga$ ORF.fasta | sed /--/d | tr "atgc" "ATGC" >> ORF_T.fasta

grep -s -B1 taa$ ORF.fasta | sed /--/d | tr "atgc" "ATGC" >> ORF_T.fasta
