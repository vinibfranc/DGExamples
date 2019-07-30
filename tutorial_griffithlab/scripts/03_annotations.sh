#!/bin/bash
set -e
set -u
set -o pipefail

echo $RNA_REFS_DIR
cd $RNA_REFS_DIR
wget http://genomedata.org/rnaseq-tutorial/annotations/GRCh38/chr22_with_ERCC92.gtf

echo $RNA_REF_GTF
less -p start_codon -S $RNA_REF_GTF

perl -ne 'if ($_ =~ /(gene_id\s\"ENSG\w+\")/){print "$1\n"}' $RNA_REF_GTF | sort | uniq | wc -l

grep ENST00000342247 $RNA_REF_GTF | less -p "exon\s" -S
