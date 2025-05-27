#$NEW 2025.n  faire 2026.0 -> avec le dernier de 2025
NEW=$1
n=$(echo $NEW | cut -f2 -d".")
AAAA=$(echo $NEW | cut -f1 -d".")
VIRUS=$2

grep -A1 -f ../src/acc2.csv $(echo $VIRUS)_ALL_$NEW.fasta | sed /--/d > $(echo $VIRUS)_$AAAA.new.fasta

bash mkNewRef90.sh $(echo $VIRUS)_2025.new.fasta $(echo $VIRUS)_$AAAA.$(( $n -1 )).fasta

cdhit-est -c 0.91 -d 0 -i toCluster.fa -o $(echo $VIRUS)Cluster.fa

cat $(echo $VIRUS)Cluster.fa $(echo $VIRUS)_$AAAA.$(( $n -1 )).fasta > $(echo $VIRUS)_$NEW.fa

cdhit-est -c 0.92 -d 0 -i $(echo $VIRUS)_$NEW.fa -o $(echo $VIRUS)_$NEW.fasta

rm toCluster.fa $(echo $VIRUS)Cluster.fa
