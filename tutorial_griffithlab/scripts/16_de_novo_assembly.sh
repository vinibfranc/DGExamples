#!/bin/bash
set -e
set -u
set -o pipefail

cd $RNA_HOME/
mkdir -p expression/stringtie/de_novo/
cd expression/stringtie/de_novo/

stringtie -p 4 -l HBR_Rep1 -o HBR_Rep1/transcripts.gtf $RNA_ALIGN_DIR/HBR_Rep1.bam
stringtie -p 4 -l HBR_Rep2 -o HBR_Rep2/transcripts.gtf $RNA_ALIGN_DIR/HBR_Rep2.bam
stringtie -p 4 -l HBR_Rep3 -o HBR_Rep3/transcripts.gtf $RNA_ALIGN_DIR/HBR_Rep3.bam

stringtie -p 4 -l UHR_Rep1 -o UHR_Rep1/transcripts.gtf $RNA_ALIGN_DIR/UHR_Rep1.bam
stringtie -p 4 -l UHR_Rep2 -o UHR_Rep2/transcripts.gtf $RNA_ALIGN_DIR/UHR_Rep2.bam
stringtie -p 4 -l UHR_Rep3 -o UHR_Rep3/transcripts.gtf $RNA_ALIGN_DIR/UHR_Rep3.bam
