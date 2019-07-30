#!/bin/bash
set -e
set -u
set -o pipefail

cd $RNA_HOME
hisat2_extract_splice_sites.py $RNA_REF_GTF > $RNA_REFS_DIR/splicesites.tsv
hisat2_extract_exons.py $RNA_REF_GTF > $RNA_REFS_DIR/exons.tsv
hisat2-build -p 8 --ss $RNA_REFS_DIR/splicesites.tsv --exon $RNA_REFS_DIR/exons.tsv $RNA_REF_FASTA $RNA_REF_INDEX


