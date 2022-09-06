#! /bin/bash

if [ -z "$*" ]; then
  echo "Arguments are missing"
  exit 0
fi

if [ -z "$2" ]; then
db="db_GROUP"
else
db=$2
fi
#Have to . is tagged (in the header name eg. "_TTV1") and if yes add sequences to the corresponding db and directory

for HEADER in $( grep ">" $1); do
  echo $HEADER | tr "_" "\n" > HEADER_TEMP   #/!\ header with "_"
  count=$(cat HEADER_TEMP | sort | uniq | wc -l)
  if [ $count -eq 2 ] ; then
    gp=$(echo $HEADER | cut -d"_" -f2)
    grep -A1 $(echo $HEADER) $1 |cut -d"_" -f1 >> $gp.fa
  elif [ $count -eq 1 ] ; then
    grep -A1 $(echo $HEADER) $1 >> no_gp.fa
  else
    grep -A1 $(echo $HEADER) $1 >> notOK_gp.fa
  fi
done