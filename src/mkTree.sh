

> ../file/TTV_TREE.fa

> ../file/TTV_B1_TREE.fa
> ../file/TTV_B2_TREE.fa
> ../file/TTV_B3_TREE.fa
> ../file/TTV_B4_TREE.fa
> ../file/TTV_B5_TREE.fa
> ../file/TTV_B6_TREE.fa
> ../file/TTV_B7_TREE.fa
> ../file/TTV_B8_TREE.fa

for i in $(ls ../db/*/*/*fasta); do echo -n ">">> ../file/TTV_TREE.fa; echo $i | cut -d "/" -f5 | cut -d "." -f1 >> ../file/TTV_TREE.fa ;grep -v ">" $i >> ../file/TTV_TREE.fa ; done

for i in $(ls ../db/B1/*/*fasta); do echo -n ">">> ../file/TTV_B1_TREE.fa; echo $i | cut -d "/" -f5 | cut -d "." -f1 >> ../file/TTV_B1_TREE.fa ;grep -v ">" $i >> ../file/TTV_B1_TREE.fa ; done

for i in $(ls ../db/B2/*/*fasta); do echo -n ">">> ../file/TTV_B2_TREE.fa; echo $i | cut -d "/" -f5 | cut -d "." -f1 >> ../file/TTV_B2_TREE.fa ;grep -v ">" $i >> ../file/TTV_B2_TREE.fa ; done

for i in $(ls ../db/B3/*/*fasta); do echo -n ">">> ../file/TTV_B3_TREE.fa; echo $i | cut -d "/" -f5 | cut -d "." -f1 >> ../file/TTV_B3_TREE.fa ;grep -v ">" $i >> ../file/TTV_B3_TREE.fa ; done

for i in $(ls ../db/B4/*/*fasta); do echo -n ">">> ../file/TTV_B4_TREE.fa; echo $i | cut -d "/" -f5 | cut -d "." -f1 >> ../file/TTV_B4_TREE.fa ;grep -v ">" $i >> ../file/TTV_B4_TREE.fa ; done

for i in $(ls ../db/B5/*/*fasta); do echo -n ">">> ../file/TTV_B5_TREE.fa; echo $i | cut -d "/" -f5 | cut -d "." -f1 >> ../file/TTV_B5_TREE.fa ;grep -v ">" $i >> ../file/TTV_B5_TREE.fa ; done

for i in $(ls ../db/B6/*/*fasta); do echo -n ">">> ../file/TTV_B6_TREE.fa; echo $i | cut -d "/" -f5 | cut -d "." -f1 >> ../file/TTV_B6_TREE.fa ;grep -v ">" $i >> ../file/TTV_B6_TREE.fa ; done

for i in $(ls ../db/B7/*/*fasta); do echo -n ">">> ../file/TTV_B7_TREE.fa; echo $i | cut -d "/" -f5 | cut -d "." -f1 >> ../file/TTV_B7_TREE.fa ;grep -v ">" $i >> ../file/TTV_B7_TREE.fa ; done

for i in $(ls ../db/B8/*/*fasta); do echo -n ">">> ../file/TTV_B8_TREE.fa; echo $i | cut -d "/" -f5 | cut -d "." -f1 >> ../file/TTV_B8_TREE.fa ;grep -v ">" $i >> ../file/TTV_B8_TREE.fa ; done



