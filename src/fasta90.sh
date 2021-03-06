#! /bin/bash
OUTPUT=sp.fa

cat $1 | tr "atgc" "ATGC" > $OUTPUT

for HEADER in $( grep ">" $OUTPUT); do
    grep -A1 $(echo $HEADER) $OUTPUT > $(echo $HEADER | cut -d ">" -f 2).fa
    echo $(echo $HEADER | cut -d ">" -f 2)" :"

       for i in $(ls db/*/*.fasta); do

          cat $(echo $i  $(echo $HEADER | cut -d ">" -f 2).fa) > file_bash_alignement_inProgress.fas
          muscle -in file_bash_alignement_inProgress.fas -clw -out file_bash_alignement.fas 2> max_len

          SAME=$(cat file_bash_alignement.fas | grep -o "*" | wc -l)

          LEN=$(cat max_len | grep "file_bash_alignement_inProgress" | cut -d " " -f6 | cut -d "," -f1)

          PTC_IDENT=$(echo $(calc $SAME/$LEN*100)| cut -d "." -f 1 | cut -d "~" -f 2)

          if [ $PTC_IDENT -gt 90 ]
          then
            echo -e "\033[1;31m" $PTC_IDENT "% identity with" $(echo $i | cut -d "." -f 1) "\033[0m"

            sed -i s/$(echo $HEADER | cut -d ">" -f 2)/$(echo $HEADER | cut -d ">" -f 2)"_"$(echo $i | cut -d "." -f 1| cut -d "/" -f 3 )/ $OUTPUT 2>/dev/null
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


