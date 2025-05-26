for i in $(ls | grep "db_" ); do cat $i/*/*.fasta > check/$i.fa; done
