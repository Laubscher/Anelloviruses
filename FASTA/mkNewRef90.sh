#usage $1 newdata 20nn.n.fasta $2-> refcluster90 20nn.n-1.fasta

sed -i s/torquevirus-hominid/torquevirus_hominid/g $1

grep ">" $1 | cut -f4 -d "," | cut -f 2 -d "=" | sort | uniq > splist.txt 


for i in $(cat splist.txt); do grep -A1 $i"," $1 | sed s/--// > all.fasta; grep -A1 $i"," $2 | sed s/--// > ref.fasta; python3 test.py >> toCluster.fa  ; done

rm splist.txt all.fasta ref.fasta
