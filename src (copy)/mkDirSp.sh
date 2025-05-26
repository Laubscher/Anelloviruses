#! /bin/bash
for i in $(grep ">" seq.fa | cut -f2 -d">") ; do mkdir $i ; grep -A1 $i seq.fa > $i/$i.fasta; done
