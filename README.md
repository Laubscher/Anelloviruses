# Anelloviruses

This repository contains scripts:
-  to create and maintain a database of complete Alphatorquevirus (TTV) ORF1 sequence annotated at species and genotype levels for metagenomic usage.
-  to analyze High-Throughput Sequencing data with the database mentioned above.


Directories:   
FASTA/        &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*.fasta version of the database*  
CSV/          &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*.csv (comma-separated values) version of the database*  
snap-db/      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*database for mapping HTS data*  
src/          &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*scripts*  


## Database format and specificity:  
version: **Alphatorquevirus_2022.1**  

FASTA format:  
>\>GENBANK,GENUS=<>,GROUP=<>,SPECIES=<>,GENOTYPE=<>,HOST=<>  
>Sequence

CSV format :
>GENBANK,GENUS,GROUP,SPECIES,GENOTYPE,HOST,Sequence


## Dependencies:

muscle can be found here:<https://www.drive5.com/muscle/downloads.htm>.  
snap-aligner can be found here:<http://snap.cs.berkeley.edu/>.  
