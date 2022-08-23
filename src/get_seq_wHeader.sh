#! /bin/bash
for i in *; do cat $i/*.fasta | sed s/\>/\>$i"_"/g >> seq.fa ; done
