# $1 version AAAA.n

bash mkFastaFinalWoGenotype.sh Betatorquevirus-A > Betatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Betatorquevirus-B >> Betatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Betatorquevirus-C >> Betatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Betatorquevirus-D >> Betatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Betatorquevirus-E >> Betatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Betatorquevirus-F >> Betatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Betatorquevirus-G >> Betatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Betatorquevirus-H >> Betatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Betatorquevirus-J >> Betatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Betatorquevirus-K >> Betatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Betatorquevirus-NA >> Betatorquevirus_$1.fasta

sed -i s/SPECIES=Betatorquevirus-/SPECIES=Betatorquevirus_/ Betatorquevirus_$1.fasta

cat  Betatorquevirus_$1.fasta | tr "%" "/" > ../FASTA/Betatorquevirus_$1.fasta
rm Betatorquevirus_$1.fasta
