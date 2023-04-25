#! /bin/bash
for i in $(grep ">" $1 | cut -f2 -d">") ; do j=$(echo $1| cut -f1 -d"." |cut -f1 -d"6" ); mkdir $j ; grep -A1 $i $1 > $j/$i.fasta; done
