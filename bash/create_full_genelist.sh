#!/usr/bin/env bash

# Define input and output files
IN=$1 # Genbank file
OUT=$2 # Tab-delimited file with features

# Header
echo -e "Type\tStart\tEnd\tLocus_tag\tOld_locus_tag\tProduct\tStrand" > $OUT

# Parse GenBank file
cat $IN | \
grep -P "(^ {5}[a-z,A-Z]+)|/old_locus_tag|/locus_tag|/pseudo|/product" | \
tr "\n" "&" | sed -e 's/ \{9\}\+/\t/g' | sed -e 's/ \{5\}/\n/g' | \
grep -P "(^CDS\t)|(^ncRNA\t)|(^rRNA\t)|(^tmRNA\t)|(^tRNA\t)" | \
sed -e '/\/old_locus_tag/! s/\/locus_tag="[a-z,A-Z,0-9,_]*"\&/&\t&/' | \
sed -e 's/ *\t */\t/g' | sed -e '/\t\/pseudo\&\t/ s/^/pseudo/' | \
sed -e 's/\t\/pseudo\&\t/\t/' | sed -e 's/)*\&\t\/locus_tag=/\t/g' \
-e 's/\/old_locus_tag=//g' -e 's/\/product=//g' | tr -d '&"' | \
sed -e 's/$/\t+/' | sed -e '/complement/ s/\t+$/\t-/' | sed -e 's/\.\./\t/g' | \
sed -e 's/complement(//' | tr -d "><" >> $OUT
