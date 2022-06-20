#! /bin/bash
#orf frame info deleted should not be here
sed -i "s/_ORF1_frAme1_/_/g" output.fa
sed -i "s/_ORF1_frAme2_/_/g" output.fa
sed -i "s/_ORF1_frAme3_/_/g" output.fa

#TODO reverse too

#Have to . is tagged (in the header name eg. "_TTV1") and if yes add sequences to the corresponding db and directory

for HEADER in $( grep ">" $1); do

  echo $HEADER | tr "_" " " > HEADER_TEMP
  count=$(cat HEADER_TEMP | sort | uniq | wc -l)
  if [ $count -eq 2 ] ; then
  grep -A1 $(echo $HEADER) $1 > OK_sp.fa
  else ;
  grep -A1 $(echo $HEADER) $1 > notOK_sp.fa
  fi
done