# $1 version AAAA.n

bash mkFastaFinalWoGenotype.sh Gammatorquevirus-A > ../FASTA/Gammatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Gammatorquevirus-B >> ../FASTA/Gammatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Gammatorquevirus-C >> ../FASTA/Gammatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Gammatorquevirus-D >> ../FASTA/Gammatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Gammatorquevirus-E >> ../FASTA/Gammatorquevirus_$1.fasta
bash mkFastaFinalWoGenotype.sh Gammatorquevirus-NA >> ../FASTA/Gammatorquevirus_$1.fasta

sed -i s/SPECIES=Gammatorquevirus-homidi/SPECIES=Gammatorquevirus_homidi/ ../FASTA/Gammatorquevirus_$1.fasta
