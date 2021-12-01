mkdir ORF1-ttv
for fasta in $(ls *.fasta);

  do sed -i s"/ //"g $fasta;

done;

for fasta in $(ls *.fasta);
 
  do wr=FALSE;
     echo ">"$fasta | cut -d "." -f 1 > ORF1-ttv/$fasta
     for line in $(cat $fasta);
        
        do  

             if [ ${line:0:1} = ">" ]; then CDS=$(echo $line | cut -d "[" -f2 | cut -d "]" -f1); echo $CDS
                                          if [ $CDS = "protein=ORF1" ]; then wr=TRUE; 
                                          elif [ $CDS = "gene=ORF1" ]; then wr=TRUE; 
                                          elif [ $CDS = "protein=ORF1a" ]; then wr=TRUE; 
                                          elif [ $CDS = "gene=ORF1a" ]; then wr=TRUE;
                                          elif [ $CDS = "protein=orf1" ]; then wr=TRUE; 
                                          elif [ $CDS = "gene=orf1" ]; then wr=TRUE;
                                          elif [ $CDS = "protein=orf1a" ]; then wr=TRUE; 
                                          elif [ $CDS = "gene=orf1a" ]; then wr=TRUE;
                                          elif [ $CDS = "protein=ORF1protein" ]; then wr=TRUE; 
                                          elif [ $CDS = "protein=ORF1aprotein" ]; then wr=TRUE;
                                          elif [ $CDS = "protein=Orf1" ]; then wr=TRUE;
                                          elif [ $CDS = "protein=orf1a" ]; then wr=TRUE;
                                          elif [ $CDS = "gene=Orf1" ]; then wr=TRUE;
                                          elif [ $CDS = "gene=Orf1a" ]; then wr=TRUE;
                                          elif [ $CDS = "protein=pORF1" ]; then wr=TRUE;
                                          else 
                                            wr=FALSE; 
                                          fi
             elif [ ${line:0:1} != ">" ]; then if [ $wr = "TRUE" ]; then echo -n $line >> ORF1-ttv/$fasta ; fi ; fi;
             
       done;

  done

for i in $(ls ORF1-ttv/*.fasta); do echo "" >> $i ; done ;
