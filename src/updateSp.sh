#! /bin/bash
for ttv in $(ls ../db/) ; do cp AnelloSP.sh ../db/$ttv/ ; cd ../db/$ttv/; bash AnelloSP.sh ../../src/$1; cd ../../src/;  done
