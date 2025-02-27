# $1 version AAAA.n

bash mkFastaFinalWoGenotype.sh Alphatorquevirus-A > Alphatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Alphatorquevirus-B >> Alphatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Alphatorquevirus-C >> Alphatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Alphatorquevirus-D >> Alphatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Alphatorquevirus-E >> Alphatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Alphatorquevirus-F >> Alphatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Alphatorquevirus-G >> Alphatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Alphatorquevirus-H >> Alphatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Alphatorquevirus-J >> Alphatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Alphatorquevirus-K >> Alphatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Alphatorquevirus-NA >> Alphatorquevirus_$1.fasta

sed -i s/SPECIES=Alphatorquevirus-/SPECIES=Alphatorquevirus_/ Alphatorquevirus_$1.fasta

cat  Alphatorquevirus_$1.fasta | tr "%" "/" > ../FASTA/Alphatorquevirus_$1.fasta
rm Alphatorquevirus_$1.fasta
