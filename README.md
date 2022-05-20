# Anelloviruses

This repository contains scripts:
-  to create and maintain a database of complete Alphatorquevirus (TTV) ORF1 sequence annotated at species and genotype levels for metagenomic usage.
-  to analyze High-Throughput Sequencing data with the database mentioned above.

Database version: Alphatorquevirus_2022.1

Database Format and specificity: 

*>FASTA/Alphatorquevirus_2022.1.fasta*

Fasta format
>\>GENBANK,GENUS=<>,GROUP=<>,SPECIES=<>,GENOTYPE=<>,HOST=<>
>Sequence

*csv/Alphatorquevirus_2022.1.csv*

csv format:
>GENBANK,GENUS,GROUP,SPECIES,GENOTYPE,HOST,Sequence

*snap-db/Alphatorquevirus_2022.1/*

index for snap-aligner

Dependencies:

muscle can be found here:<https://www.drive5.com/muscle/downloads.htm>.  
snap-aligner can be found here:<http://snap.cs.berkeley.edu/>.  
Samtools can be found here:<http://www.htslib.org/download/>.  
Bedtools can be found here:<https://bedtools.readthedocs.io/en/latest/>.  
tagdust can be found here:<http://tagdust.sourceforge.net/>.  
trimmomatic can be found here:<http://www.usadellab.org/cms/?page=trimmomatic>.  
