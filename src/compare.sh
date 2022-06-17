#! /bin/bash

if [ -z "$*" ]; then
  echo "Arguments are missing (compare.sh -h for help)"
  exit 0
fi

if [ -z "$2" ]; then
PCT=69                # 69% for sp. determination Taxonomic update for mammalian anelloviruses (family Anelloviridae)
else                  # Archives of Virology (2021) 166:2943–2953 https://doi.org/10.1007/s00705-021-05192-x
PCT=$2
fi

if [ -z "$3" ]; then
db="db_Alphatorquevirus"
else
db=$3
fi

while getopts ":h" opt; do
	case $opt in
                 h)
                       echo "Usage : compare.sh -options(-h display this help) [fastafile.fa] [% nucl identity] [db path]"
                       exit 0
                       ;;

		\?)
			echo "Unknown option $OPTARG" >&2
			;;
	esac
done




OUTPUT=output.fa
#input a fasta file
#
#output header has _tagged dir of matching sequence

cat $1 | tr "atgc" "ATGC" > $OUTPUT

for HEADER in $( grep ">" $OUTPUT); do
    grep -A1 $(echo $HEADER) $OUTPUT > $(echo $HEADER | cut -d ">" -f 2).fa
    echo $(echo $HEADER | cut -d ">" -f 2)" :"

       for i in $(ls $db/*/*.fasta); do
          cat $(echo $i  $(echo $HEADER | cut -d ">" -f 2).fa) > file_bash_alignement_inProgress.fas
          muscle -in file_bash_alignement_inProgress.fas -clw -out file_bash_alignement.fas 2> max_len

          SAME=$(cat file_bash_alignement.fas | grep -o "*" | wc -l)

          LEN=$(cat max_len | grep "file_bash_alignement_inProgress" | cut -d " " -f6 | cut -d "," -f1)

          PTC_IDENT=$(echo $(calc $SAME/$LEN*100)| cut -d "." -f 1 | cut -d "~" -f 2)

          if [ $PTC_IDENT -gt $PCT ]
          then
            echo -e "\033[1;31m" $PTC_IDENT "% identity with" $(echo $i | cut -d "." -f 1) "\033[0m"

            sed -i s/$(echo $HEADER | cut -d ">" -f 2)/$(echo $HEADER | cut -d ">" -f 2)"_"$(echo $i | cut -d "." -f 1| cut -d "/" -f 2 )/ $OUTPUT 2>/dev/null
          else
            echo " "$PTC_IDENT "% identity with" $(echo $i | cut -d "." -f 1)
          fi
       done
    echo

    rm file_bash_alignement_inProgress.fas
    rm $(echo $HEADER | cut -d ">" -f 2).fa
    rm max_len
    rm file_bash_alignement.fas
done


