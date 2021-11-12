> TTV_TREE.fa
for i in $(ls TTV*/*fasta); do echo -n ">">> TTV_TREE.fa; echo $i | cut -d "/" -f2 | cut -d "." -f1 >> TTV_TREE.fa ;grep -v ">" $i >> TTV_TREE.fa ; done
