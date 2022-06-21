#! /bin/bash

if [ -z "$*" ]; then
  echo "Arguments are missing"
  exit 0
fi

if [ -z "$2" ]; then
db="db_Alphatorquevirus"
else
db=$2
fi

#Have to . is tagged (in the header name eg. "_TTV1") and if yes add sequences to the corresponding db and directory

for HEADER in $( grep ">" $1); do
  echo $HEADER | tr "_" "\n" > HEADER_TEMP   #/!\ header with "_"
  count=$(cat HEADER_TEMP | sort | uniq | wc -l)
  if [ $count -eq 2 ] ; then
    sp=$(echo $HEADER | cut -d"_" -f2)
    ac=$(echo $HEADER | cut -d"_" -f1 | cut -d">" -f2 )
    grep -A1 $(echo $HEADER) $1 |cut -d"_" -f1 > $db/$sp/$ac.fasta
  elif [ $count -eq 1 ] ; then
    grep -A1 $(echo $HEADER) $1 >> no_sp.fa
  else
    grep -A1 $(echo $HEADER) $1 >> notOK_sp.fa
  fi
done