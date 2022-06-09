#!/bin/bash

NAME=$(echo $1 |rev |cut -d "/" -f 1 | rev) 
OUTPUT=$(echo $4)
TEMPNAME=$(echo $OUTPUT/$NAME)
CPU=$(echo $2)
VIRUSDB=$(echo $3)

TIME=$(date +%s)

cp $1".R1.fastq.gz" $TEMPNAME"_R1_001.fastq.gz" 
cp $1".R2.fastq.gz" $TEMPNAME"_R2_001.fastq.gz"

echo "Copying :"$(($(date +%s) - $TIME))" sec."
TIME=$(date +%s)

gunzip $(echo $TEMPNAME)_R1_001.fastq.gz
gunzip $(echo $TEMPNAME)_R2_001.fastq.gz

echo "Ungziping :"$(($(date +%s) - $TIME))" sec."
TIME=$(date +%s)


FQ1=$(echo $TEMPNAME)_R1_001.fastq

FQ2=$(echo $TEMPNAME)_R2_001.fastq


tagdust -t $CPU -dust 16 -1 R:N $FQ1 $FQ2 -o $TEMPNAME

rm $(echo $TEMPNAME)_logfile.txt $(echo $TEMPNAME)_un_READ1.fq $(echo $TEMPNAME)_un_READ2.fq
  
echo "Low-Complexity Filtering :"$(($(date +%s) - $TIME))" sec."
TIME=$(date +%s)

mv $(echo $TEMPNAME)_READ1.fq $(echo $TEMPNAME)_R1_001.fastq

mv $(echo $TEMPNAME)_READ2.fq $(echo $TEMPNAME)_R2_001.fastq

TIME=$(date +%s)

bash virus_mapping $TEMPNAME $VIRUSDB $CPU 

echo "Virus Mapping :"$(($(date +%s) - $TIME))" sec."
TIME=$(date +%s)

# /!\ >= e+06 value in .bga have to be convert to numerical

python3 metrics.py $TEMPNAME $(echo $3)/ref.txt files/black_list.txt files/colors.txt

echo "Metrics Calculating :"$(($(date +%s) - $TIME))" sec."

rm $(echo $TEMPNAME)*.bam # may want to keep it

exit 0
