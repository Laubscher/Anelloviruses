# Anelloviruses

This repository contains scripts:
-  to create and maintain a database of complete Alphatorquevirus (TTV) ORF1 sequence annotated at species and genotype levels for metagenomic usage.
-  to analyze High-Throughput Sequencing data with the database mentioned above.


Directories: 
FASTA/      *.fasta version of the database*  
CSV/        *.csv (comma-separated values) version of the database*  
snap-db/    *database for mapping HTS data*  
src/        *scripts*  


Database Format and specificity: 
version: Alphatorquevirus_2022.1

FASTA format:
>\>GENBANK,GENUS=<>,GROUP=<>,SPECIES=<>,GENOTYPE=<>,HOST=<>  
>Sequence

CSV format :
>GENBANK,GENUS,GROUP,SPECIES,GENOTYPE,HOST,Sequence

Dependencies:

muscle can be found here:<https://www.drive5.com/muscle/downloads.htm>.  
snap-aligner can be found here:<http://snap.cs.berkeley.edu/>.  
