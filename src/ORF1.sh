#! /bin/bash

FILE=$(echo $1) ;

   cat $FILE | tr "\n" "!" > $FILE.temp               
   sed -i "s/>/\n/g" $FILE.temp                       
   sed -i "s/A/A_/g" $FILE.temp
   sed -i "s/T/T_/g" $FILE.temp
   sed -i "s/G/G_/g" $FILE.temp
   sed -i "s/C/C_/g" $FILE.temp
   sed -i "s/N/N_/g" $FILE.temp 2>/dev/null
   
   for ENTRY in $( cat $FILE.temp) ;
     do 
     echo $(echo $ENTRY | cut -d "!" -f 1 | cut -d "|" -f 1)
     rm $FILE.FRAME1.temp $FILE.FRAME2.temp $FILE.FRAME3.temp $FILE.FRAME1.rev.temp $FILE.FRAME2.rev.temp $FILE.FRAME3.rev.temp 2>/dev/null
     SEQ=$( echo $ENTRY | cut -d "!" -f 2 )
     I=0     

     for NUCL in $( echo $SEQ | tr "_" " ") ;           #For each nucleotides
       do
       I=$(( $I + 1))
       echo -n $NUCL >> $FILE.FRAME1.temp
       echo -n $NUCL >> $FILE.FRAME2.temp
       echo -n $NUCL >> $FILE.FRAME3.temp   

       if [ $(($I%3)) -eq 0 ]
         then 
         echo -n " " >> $FILE.FRAME1.temp
       fi
       if [ $(($I%3)) -eq 2 ]
         then 
         echo -n " " >> $FILE.FRAME2.temp
       fi
       if [ $(($I%3)) -eq 1 ]
         then 
         echo -n " " >> $FILE.FRAME3.temp
       fi
     done

     I=0     

     SEQ=$( echo $ENTRY | cut -d "!" -f 2 | tr "ATGC"  "atgc" | tr "atgc" "TACG" | rev )

     for NUCL in $( echo $SEQ | tr "_" " ") ;
       do
       I=$(( $I + 1))
       echo -n $NUCL >> $FILE.FRAME1.rev.temp
       echo -n $NUCL >> $FILE.FRAME2.rev.temp
       echo -n $NUCL >> $FILE.FRAME3.rev.temp


       if [ $(($I%3)) -eq 0 ]
         then 
         echo -n " " >> $FILE.FRAME1.rev.temp
       fi
       if [ $(($I%3)) -eq 2 ]
         then 
         echo -n " " >> $FILE.FRAME2.rev.temp
       fi
       if [ $(($I%3)) -eq 1 ]
         then 
         echo -n " " >> $FILE.FRAME3.rev.temp
       fi
     
    done

    sed -i "s/ATG/xatg/g" $FILE.FRAME1.temp 2>/dev/null
    sed -i "s/ACG/xacg/g" $FILE.FRAME1.temp 2>/dev/null    # for TTV with ACG START codon aka TTV 27 ,28 ,29
    sed -i "s/TGA/tgaX/g" $FILE.FRAME1.temp 2>/dev/null  
    sed -i "s/TAG/tagX/g" $FILE.FRAME1.temp 2>/dev/null 
    sed -i "s/TAA/taaX/g" $FILE.FRAME1.temp 2>/dev/null
    sed -i "s/ //g" $FILE.FRAME1.temp 2>/dev/null

    sed -i "s/ATG/xatg/g" $FILE.FRAME2.temp 2>/dev/null
    sed -i "s/ACG/xacg/g" $FILE.FRAME2.temp 2>/dev/null
    sed -i "s/TGA/tgaX/g" $FILE.FRAME2.temp 2>/dev/null
    sed -i "s/TAG/tagX/g" $FILE.FRAME2.temp 2>/dev/null
    sed -i "s/TAA/taaX/g" $FILE.FRAME2.temp 2>/dev/null
    sed -i "s/ //g" $FILE.FRAME2.temp 2>/dev/null

    sed -i "s/ATG/xatg/g" $FILE.FRAME3.temp 2>/dev/null
    sed -i "s/ACG/xacg/g" $FILE.FRAME3.temp 2>/dev/null 
    sed -i "s/TGA/tgaX/g" $FILE.FRAME3.temp 2>/dev/null
    sed -i "s/TAG/tagX/g" $FILE.FRAME3.temp 2>/dev/null
    sed -i "s/TAA/taaX/g" $FILE.FRAME3.temp 2>/dev/null
    sed -i "s/ //g" $FILE.FRAME3.temp 2>/dev/null

    sed -i "s/ATG/xatg/g" $FILE.FRAME1.rev.temp 2>/dev/null
    sed -i "s/ACG/xacg/g" $FILE.FRAME1.rev.temp 2>/dev/null
    sed -i "s/TGA/tgaX/g" $FILE.FRAME1.rev.temp 2>/dev/null
    sed -i "s/TAG/tagX/g" $FILE.FRAME1.rev.temp 2>/dev/null
    sed -i "s/TAA/taaX/g" $FILE.FRAME1.rev.temp 2>/dev/null
    sed -i "s/ //g" $FILE.FRAME1.rev.temp 2>/dev/null

    sed -i "s/ATG/xatg/g" $FILE.FRAME2.rev.temp 2>/dev/null
    sed -i "s/ACG/xacg/g" $FILE.FRAME2.rev.temp 2>/dev/null
    sed -i "s/TGA/tgaX/g" $FILE.FRAME2.rev.temp 2>/dev/null
    sed -i "s/TAG/tagX/g" $FILE.FRAME2.rev.temp 2>/dev/null
    sed -i "s/TAA/taaX/g" $FILE.FRAME2.rev.temp 2>/dev/null
    sed -i "s/ //g" $FILE.FRAME2.rev.temp 2>/dev/null

    sed -i "s/ATG/xatg/g" $FILE.FRAME3.rev.temp 2>/dev/null
    sed -i "s/ACG/xacg/g" $FILE.FRAME3.rev.temp 2>/dev/null
    sed -i "s/TGA/tgaX/g" $FILE.FRAME3.rev.temp 2>/dev/null
    sed -i "s/TAG/tagX/g" $FILE.FRAME3.rev.temp 2>/dev/null
    sed -i "s/TAA/taaX/g" $FILE.FRAME3.rev.temp 2>/dev/null
    sed -i "s/ //g" $FILE.FRAME3.rev.temp 2>/dev/null

  
    for i in $( cat $FILE.FRAME1.temp | tr "X" " ");                    # cut at all stop condon
      do
        
        i=$(echo $i | sed s/$(echo $i |cut -d "x" -f 1)//  2>/dev/null)           # all what is before the firts start codon is erase
        i=$(echo $i | sed s/x//g 2>/dev/null)                                     # clear all the markt for start codon      
        
        #ri=$(echo $i | rev)

        #if [-z $(echo ${ri:0:3} | sed s/tga// | sed | s/tga// | sed s/tga//) ]

        if [ ${#i} -gt 1799 ]
          then
          echo ">"$(echo $ENTRY | cut -d "!" -f 1 | cut -d "|" -f 1)"_ORF1_frame1" >> ORF.fasta
          echo $i  >> ORF.fasta
        fi
    done   

    for i in $( cat $FILE.FRAME2.temp | tr "X" " ");
      do
        
        i=$(echo $i | sed s/$(echo $i |cut -d "x" -f 1)//  2>/dev/null) 
        i=$(echo $i | sed s/x//g 2>/dev/null) 

        if [ ${#i} -gt 1799 ]
          then
          echo ">"$(echo $ENTRY | cut -d "!" -f 1 | cut -d "|" -f 1 )"_ORF1_frame2" >> ORF.fasta
          echo $i  >> ORF.fasta
        fi
    done

    for i in $( cat $FILE.FRAME3.temp | tr "X" " ");
      do
        
        i=$(echo $i | sed s/$(echo $i |cut -d "x" -f 1)//  2>/dev/null) 
        i=$(echo $i | sed s/x//g 2>/dev/null) 

        if [ ${#i} -gt 1799 ]
          then
          echo ">"$(echo $ENTRY | cut -d "!" -f 1 | cut -d "|" -f 1)"_ORF1_frame3" >> ORF.fasta
          echo $i  >> ORF.fasta
        fi
    done   

    for i in $( cat $FILE.FRAME1.rev.temp | tr "X" " ");
      do
        
        i=$(echo $i | sed s/$(echo $i |cut -d "x" -f 1)//  2>/dev/null) 
        i=$(echo $i | sed s/x//g 2>/dev/null) 

        if [ ${#i} -gt 1799 ]
          then
          echo ">"$(echo $ENTRY | cut -d "!" -f 1 | cut -d "|" -f 1)"_ORF1_frame1_rev" >> ORF.fasta
          echo $i  >> ORF.fasta
        fi
    done   

    for i in $( cat $FILE.FRAME2.rev.temp | tr "X" " ");
      do
        
        i=$(echo $i | sed s/$(echo $i |cut -d "x" -f 1)//  2>/dev/null) 
        i=$(echo $i | sed s/x//g 2>/dev/null) 

        if [ ${#i} -gt 1799 ]
          then
          echo ">"$(echo $ENTRY | cut -d "!" -f 1 | cut -d "|" -f 1)"_ORF1_frame2_rev" >> ORF.fasta
          echo $i  >> ORF.fasta
        fi
    done

    for i in $( cat $FILE.FRAME3.rev.temp | tr "X" " ");
      do
        
        i=$(echo $i | sed s/$(echo $i |cut -d "x" -f 1)//  2>/dev/null) 
        i=$(echo $i | sed s/x//g 2>/dev/null) 

        if [ ${#i} -gt 1799 ]
          then
          echo ">"$(echo $ENTRY | cut -d "!" -f 1 | cut -d "|" -f 1)"_ORF1_frame3_rev" >> ORF.fasta
          echo $i  >> ORF.fasta
        fi
    done 

  done
  rm $FILE.temp $FILE.FRAME1.temp $FILE.FRAME2.temp $FILE.FRAME3.temp
