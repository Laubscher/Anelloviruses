> TTV_TREE.fa
for i in $(ls ../db/*/*fasta); do echo -n ">">> ../file/TTV_TREE.fa; echo $i | cut -d "/" -f4 | cut -d "." -f1 >> ../file/TTV_TREE.fa ;grep -v ">" $i >> ../file/TTV_TREE.fa ; done
