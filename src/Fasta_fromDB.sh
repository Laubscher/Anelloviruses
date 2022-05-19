#Get FASTA OF ORF1 TTV 
#$1 = database.csv

sed '$d' $1 > temp_fasta.

# Enlever le header et remplacer les espaces par _

sed -i '1d' temp_fasta.
sed -i s"/ /_/"g temp_fasta.

#Recrer l'header et la sequence

prefix=$(echo $1 | cut -d"." -f1)
>$prefix.fa
for line in $(cat temp_fasta.);
  do Header1=">"$(echo $line | cut -f1 -d",")"," ;  
     Header2="GENUS="$(echo $line | cut -f2 -d",")"," ;
     Header3="GROUP="$(echo $line | cut -f3 -d",")"," ;
     Header4="SPECIES="$(echo $line | cut -f4 -d",")"," ;
     Header5="GENOTYPE="$(echo $line | cut -f5 -d",")"," ;
     Header6="HOST="$(echo $line | cut -f6 -d",") ;
     seq=$(echo $line | cut -f7 -d",") ;
     echo $Header1 $Header2 $Header3 $Header4 $Header5 $Header6 >> $prefix.fa
     echo $seq >> $prefix.fa
                                                  done ;
sed -i s"/ //"g $prefix.fa
rm temp_fasta.			

