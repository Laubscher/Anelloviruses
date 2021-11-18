#! /bin/bash
#$1 fasta.fa query to assigne

outputfile=not_assigned$(date | tr " :" "__").fa
> $outputfile
ls *.fasta| cut -f5 -d "_" | cut -f1 -d "." > listN.txt
grep -vf listN.txt $1 > newN.fatemp
for HEADER in $( grep ">" newN.fatemp); do
    grep -A1 $(echo $HEADER) newN.fatemp >> $outputfile
done
rm listN.txt
rm newN.fatemp
