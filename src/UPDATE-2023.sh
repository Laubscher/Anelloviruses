#script for updating database

#  -> sequence file -> "update_2023"

#Alphatorquevirus

python3 getFastaFromGB.py A > updating_2023         # Ask genbank with ttv taxid 

bash updateGp.sh db_Group-Alpha  # no parameters => + ORF1 finding

rm updating_2023
mv Alpha-A.fa updating_2023
bash update.sh db_Alphatorquevirus-A 68 TRUE         # get ORF1 for each sequences 
mv no_sp.fa UNCLASS/UNCLASS-Alpha-A.fa


rm updating_2023
mv Alpha-B.fa updating_2023
bash update.sh db_Alphatorquevirus-B 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Alpha-B.fa


rm updating_2023
mv Alpha-C.fa updating_2023
bash update.sh db_Alphatorquevirus-C 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Alpha-C.fa


rm updating_2023
mv Alpha-D.fa updating_2023
bash update.sh db_Alphatorquevirus-D 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Alpha-D.fa


rm updating_2023
mv Alpha-E.fa updating_2023
bash update.sh db_Alphatorquevirus-E 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Alpha-E.fa


rm updating_2023
mv Alpha-F.fa updating_2023
bash update.sh db_Alphatorquevirus-F 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Alpha-F.fa


rm updating_2023
mv Alpha-G.fa updating_2023
bash update.sh db_Alphatorquevirus-G 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Alpha-G.fa


rm updating_2023
mv Alpha-H.fa updating_2023
bash update.sh db_Alphatorquevirus-H 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Alpha-H.fa


rm updating_2023
mv Alpha-J.fa updating_2023
bash update.sh db_Alphatorquevirus-J 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Alpha-J.fa


rm updating_2023
mv Alpha-K.fa updating_2023
bash update.sh db_Alphatorquevirus-K 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Alpha-K.fa


rm updating_2023
mv Alpha-NA.fa updating_2023
bash update.sh db_Alphatorquevirus-NA 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Alpha-NA.fa

#BETA

python3 getFastaFromGB.py B > updating_2023      # Ask genbank with ttmv taxid

bash updateGp.sh db_Group-Beta

rm updating_2023
mv Beta-A.fa updating_2023
bash update.sh db_Betatorquevirus-A 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Beta-A.fa

rm updating_2023
mv Beta-B.fa updating_2023
bash update.sh db_Betatorquevirus-B 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Beta-B.fa

rm updating_2023
mv Beta-B2.fa updating_2023
bash update.sh db_Betatorquevirus-NA 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Beta-NA.fa

rm updating_2023
mv Beta-C.fa updating_2023
bash update.sh db_Betatorquevirus-C 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Beta-C.fa

rm updating_2023
mv Beta-D.fa updating_2023
bash update.sh db_Betatorquevirus-D 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Beta-D.fa

rm updating_2023
mv Beta-E.fa updating_2023
bash update.sh db_Betatorquevirus-E 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Beta-E.fa

rm updating_2023
mv Beta-F.fa updating_2023
bash update.sh db_Betatorquevirus-F 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Beta-F.fa

rm updating_2023
mv Beta-G.fa updating_2023
bash update.sh db_Betatorquevirus-G 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Beta-G.fa

rm updating_2023
mv Beta-H.fa updating_2023
bash update.sh db_Betatorquevirus-H 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Beta-H.fa

rm updating_2023
mv Beta-J.fa updating_2023
bash update.sh db_Betatorquevirus-J 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Beta-J.fa

rm updating_2023
mv Beta-K.fa updating_2023
bash update.sh db_Betatorquevirus-K 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Beta-K.fa

#GAMMA
python3 getFastaFromGB.py G > updating_2023      # Ask genbank with ttmdv taxid

bash updateGp.sh db_Group-Gamma

rm updating_2023
mv Gamma-A.fa updating_2023
bash update.sh db_Gammatorquevirus-A 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Gamma-A.fa

rm updating_2023
mv Gamma-B.fa updating_2023
bash update.sh db_Gammatorquevirus-B 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Gamma-B.fa

rm updating_2023
mv Gamma-C.fa updating_2023
bash update.sh db_Gammatorquevirus-C 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Gamma-C.fa

rm updating_2023
mv Gamma-D.fa updating_2023
bash update.sh db_Gammatorquevirus-D 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Gamma-D.fa

rm updating_2023
mv Gamma-E.fa updating_2023
bash update.sh db_Gammatorquevirus-E 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Gamma-E.fa

#OMEGA
python3 getFastaFromGB.py O > updating_2023

bash update.sh db_Omegatorquevirus
mv no_sp.fa UNCLASS/UNCLASS-Omega.fa

#HE
python3 getFastaFromGB.py H > updating_2023

bash update.sh db_Hetorquevirus
mv no_sp.fa UNCLASS/UNCLASS-He.fa

# -1 For others unclassified sequences with custums research terms (Host Human)

python3 getFastaFromGB.py -1 > updating_2023

#have to checked all genera/groups 

mv no_gp.fa no_ABG.fa

bash updateGp.sh db_Group-Alpha  # no parameters => + ORF1 finding

mv no_gp.fa updating_2023
bash updateGp.sh db_Group-Beta 63 TRUE

mv no_gp.fa updating_2023
bash updateGp.sh db_Group-Gamma 63 TRUE

mv no_gp.fa updating_2023

cat no_ABG.fa >> updating_2023
rm no_ABG.fa

bash update.sh db_Hetorquevirus 68 TRUE
mv no_sp.fa updating_2023

bash update.sh db_Omegatorquevirus 68 TRUE
mv no_sp.fa updating_2023

bash update.sh db_Lamedtorquevirus 68 TRUE
mv no_sp.fa updating_2023
bash update.sh db_Memtorquevirus 68 TRUE
mv no_sp.fa updating_2023
bash update.sh db_Samektorquevirus 68 TRUE
mv no_sp.fa updating_2023
bash update.sh db_Yodtorquevirus 68 TRUE
mv no_sp.fa updating_2023

bash update.sh db_Epsilontorquevirus 68 TRUE
mv no_sp.fa updating_2023
bash update.sh db_Zetatorquevirus 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-GROUPS.fa

####################################################################
# inside TTV TTMV TTMDV groups sp assignation 

rm updating_2023
mv Alpha-A.fa updating_2023
bash update.sh db_AlphatorquevirusA 68 TRUE        
mv no_sp.fa UNCLASS/UNCLASS-Alpha-A.fa


rm updating_2023
mv Alpha-B.fa updating_2023
bash update.sh db_AlphatorquevirusB 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Alpha-B.fa


rm updating_2023
mv Alpha-C.fa updating_2023
bash update.sh db_AlphatorquevirusC 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Alpha-C.fa


rm updating_2023
mv Alpha-D.fa updating_2023
bash update.sh db_AlphatorquevirusD 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Alpha-D.fa


rm updating_2023
mv Alpha-E.fa updating_2023
bash update.sh db_AlphatorquevirusE 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Alpha-E.fa


rm updating_2023
mv Alpha-F.fa updating_2023
bash update.sh db_AlphatorquevirusF 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Alpha-F.fa


rm updating_2023
mv Alpha-G.fa updating_2023
bash update.sh db_AlphatorquevirusG 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Alpha-G.fa


rm updating_2023
mv Alpha-H.fa updating_2023
bash update.sh db_AlphatorquevirusH 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Alpha-H.fa


rm updating_2023
mv Alpha-J.fa updating_2023
bash update.sh db_AlphatorquevirusJ 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Alpha-J.fa


rm updating_2023
mv Alpha-K.fa updating_2023
bash update.sh db_AlphatorquevirusK 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Alpha-K.fa


rm updating_2023
mv Alpha-NA.fa updating_2023
bash update.sh db_AlphatorquevirusNA 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Alpha-NA.fa


rm updating_2023
mv Beta-A.fa updating_2023
bash update.sh db_Betatorquevirus-A 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Beta-A.fa

rm updating_2023
mv Beta-B.fa updating_2023
bash update.sh db_Betatorquevirus-B 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Beta-B.fa

rm updating_2023
mv Beta-B2.fa updating_2023
bash update.sh db_Betatorquevirus-NA 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Beta-NA.fa

rm updating_2023
mv Beta-C.fa updating_2023
bash update.sh db_Betatorquevirus-C 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Beta-C.fa

rm updating_2023
mv Beta-D.fa updating_2023
bash update.sh db_Betatorquevirus-D 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Beta-D.fa

rm updating_2023
mv Beta-E.fa updating_2023
bash update.sh db_Betatorquevirus-E 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Beta-E.fa

rm updating_2023
mv Beta-F.fa updating_2023
bash update.sh db_Betatorquevirus-F 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Beta-F.fa

rm updating_2023
mv Beta-G.fa updating_2023
bash update.sh db_Betatorquevirus-G 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Beta-G.fa

rm updating_2023
mv Beta-H.fa updating_2023
bash update.sh db_Betatorquevirus-H 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Beta-H.fa

rm updating_2023
mv Beta-J.fa updating_2023
bash update.sh db_Betatorquevirus-J 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Beta-J.fa

rm updating_2023
mv Beta-K.fa updating_2023
bash update.sh db_Betatorquevirus-K 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Beta-K.fa

rm updating_2023
mv Gamma-A.fa updating_2023
bash update.sh db_GammatorquevirusA 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Gamma-A.fa

rm updating_2023
mv Gamma-B.fa updating_2023
bash update.sh db_GammatorquevirusB 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Gamma-B.fa

rm updating_2023
mv Gamma-C.fa updating_2023
bash update.sh db_GammatorquevirusC 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Gamma-C.fa

rm updating_2023
mv Gamma-D.fa updating_2023
bash update.sh db_GammatorquevirusD 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Gamma-D.fa

rm updating_2023
mv Gamma-E.fa updating_2023
bash update.sh db_GammatorquevirusE 68 TRUE
mv no_sp.fa UNCLASS/UNCLASS-Gamma-E.fa

#auto test UNCLASS/UNCLASS-GROUPS.fa for all sp?
#same forgp_not_ok.fa ?





