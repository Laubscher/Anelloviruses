#! /bin/bash
for i in $(grep ">" no_sp.fa | cut -f2 -d">") ; do mkdir $i ; grep -A1 $i no_sp.fa > $i/$i.fasta; done
