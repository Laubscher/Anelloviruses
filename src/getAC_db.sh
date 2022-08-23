#! /bin/bash
cat */*.fasta | grep ">" | cut -f2 -d">"
