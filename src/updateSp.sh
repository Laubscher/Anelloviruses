for ttv in $(ls ../db/) ; do cp AnelloSP.sh ../db/$ttv/ ; cd ../db/$ttv/; bash AnelloSP.sh ../../file/NEW_CDS.fa; cd ../../src/;  done
