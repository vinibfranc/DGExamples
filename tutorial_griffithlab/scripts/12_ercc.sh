
#!/bin/bash
set -e
set -u
set -o pipefail

cd $RNA_HOME/de/ballgown/ref_only
wget https://raw.githubusercontent.com/griffithlab/rnaseq_tutorial/master/scripts/Tutorial_ERCC_DE.R
chmod +x Tutorial_ERCC_DE.R
./Tutorial_ERCC_DE.R $RNA_HOME/expression/htseq_counts/ERCC_Controls_Analysis.txt $RNA_HOME/de/ballgown/ref_only/UHR_vs_HBR_gene_results.tsv

