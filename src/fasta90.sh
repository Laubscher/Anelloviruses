#! /bin/bash
OUTPUT=sp.fa

cat $1 | tr "atgc" "ATGC" > $OUTPUT

for HEADER in $( grep ">" $OUTPUT); do
    grep -A1 $(echo $HEADER) $OUTPUT > $(echo $HEADER | cut -d ">" -f 2).fa
    echo $(echo $HEADER | cut -d ">" -f 2)" :"

       for i in $(ls *.fasta); do

          cat $(echo $i  $(echo $HEADER | cut -d ">" -f 2).fa) > file_bash_alignement_inProgress.fas
          SAME=$(muscle -in file_bash_alignement_inProgress.fas -clw -quiet | grep -o "*" | wc -l)

          LEN=$(muscle -in file_bash_alignement_inProgress.fas -msf -quiet  | grep MSF | cut -d " " -f4)".0"

          PTC_IDENT=$(echo $(calc $SAME/$LEN*100)| cut -d "." -f 1 | cut -d "~" -f 2)
          
          if [ $PTC_IDENT -gt 90 ]  
          then
            echo -e "\033[1;31m" $PTC_IDENT "% identity with" $(echo $i | cut -d "." -f 1) "\033[0m"

            sed -i s/$(echo $HEADER | cut -d ">" -f 2)/$(echo $HEADER | cut -d ">" -f 2)"_"$(echo $i | cut -d "." -f 1)/ $OUTPUT 2>/dev/null
          else
            echo " "$PTC_IDENT "% identity with" $(echo $i | cut -d "." -f 1)
          fi
       done
    echo
    
    rm  file_bash_alignement_inProgress.fas
    rm  $(echo $HEADER | cut -d ">" -f 2).fa
done



