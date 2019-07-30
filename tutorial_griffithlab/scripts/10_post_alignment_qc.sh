#!/bin/bash
set -e
set -u
set -o pipefail

cd $RNA_ALIGN_DIR
# samtools view -H UHR.bam
# samtools view UHR.bam | head

# Try requiring that alignments are paired and mapped in a proper pair (=3). Also filter out alignments that are unmapped, the mate is unmapped, and not primary alignment (=268)
# samtools view -f 3 -F 268 UHR.bam

# samtools view -f 1024 UHR.bam

# Summary of the alignment
cd $RNA_ALIGN_DIR
samtools flagstat UHR.bam
samtools flagstat HBR.bam

# QC report using RSeQC
cd $RNA_HOME/data/
wget http://genomedata.org/rnaseq-tutorial/RSeQC.zip
unzip RSeQC.zip
cd RSeQC/
gunzip hg19_UCSC_knownGene.bed.gz

# RSeQC commands
bam_stat.py -i Pairend_nonStrandSpecific_36mer_Human_hg19.bam
clipping_profile.py -i Pairend_nonStrandSpecific_36mer_Human_hg19.bam -o tutorial -s "PE"
geneBody_coverage.py -r hg19_UCSC_knownGene.bed -i Pairend_nonStrandSpecific_36mer_Human_hg19.bam -o tutorial
infer_experiment.py -r hg19_UCSC_knownGene.bed -i Pairend_nonStrandSpecific_36mer_Human_hg19.bam
inner_distance.py -r hg19_UCSC_knownGene.bed -i Pairend_nonStrandSpecific_36mer_Human_hg19.bam -o tutorial
junction_annotation.py -r hg19_UCSC_knownGene.bed -i Pairend_nonStrandSpecific_36mer_Human_hg19.bam -o tutorial
junction_saturation.py -r hg19_UCSC_knownGene.bed -i Pairend_nonStrandSpecific_36mer_Human_hg19.bam -o tutorial
read_distribution.py -r hg19_UCSC_knownGene.bed -i Pairend_nonStrandSpecific_36mer_Human_hg19.bam
read_duplication.py -i Pairend_nonStrandSpecific_36mer_Human_hg19.bam -o tutorial
read_GC.py -i Pairend_nonStrandSpecific_36mer_Human_hg19.bam -o tutorial
read_NVC.py -i Pairend_nonStrandSpecific_36mer_Human_hg19.bam -o tutorial
read_quality.py -i Pairend_nonStrandSpecific_36mer_Human_hg19.bam -o tutorial
ls *.pdf
