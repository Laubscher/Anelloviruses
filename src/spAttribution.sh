#! /bin/bash
#orf frame info deleted should not be here
sed -i "s/_ORF1_frAme1_/_/g" output.fa
sed -i "s/_ORF1_frAme2_/_/g" output.fa
sed -i "s/_ORF1_frAme3_/_/g" output.fa

#TODO reverse too

#Have to check if exactly 1 sp. is tagged (in the header name eg. "_TTV1") and if yes add sequences to the corresponding db and directory