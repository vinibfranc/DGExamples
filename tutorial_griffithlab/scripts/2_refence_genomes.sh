#!/bin/bash
set -e
set -u
set -o pipefail

cd $RNA_HOME
echo $RNA_REFS_DIR
mkdir -p $RNA_REFS_DIR

cd $RNA_REFS_DIR
wget http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa
ls 

head chr22_with_ERCC92.fa

wc chr22_with_ERCC92.fa

head -n 425000 chr22_with_ERCC92.fa | tail

cat chr22_with_ERCC92.fa | grep -v ">" | perl -ne 'chomp $_; $bases{$_}++ for split //; if (eof){print "$_ $bases{$_}\n" for sort keys %bases}'
