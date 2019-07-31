#!/bin/bash
set -e
set -u
set -o pipefail

# In this module we will run Stringtie in two additional modes: (1) 'reference guided' mode and (2) 'de novo' mode. Stringtie can predict the transcripts present in each library with or without help from knowledge of known transcripts.

# To run Stringtie in 'reference guided' mode: use the '-G' option WITHOUT '-e'

# To run Stringtie in 'de novo' mode do NOT specify either of the '-G' OR '-e' options.

cd $RNA_HOME/
mkdir -p expression/stringtie/ref_guided/
cd expression/stringtie/ref_guided/

stringtie -p 4 -G $RNA_REF_GTF -l HBR_Rep1 -o HBR_Rep1/transcripts.gtf $RNA_ALIGN_DIR/HBR_Rep1.bam
stringtie -p 4 -G $RNA_REF_GTF -l HBR_Rep2 -o HBR_Rep2/transcripts.gtf $RNA_ALIGN_DIR/HBR_Rep2.bam
stringtie -p 4 -G $RNA_REF_GTF -l HBR_Rep3 -o HBR_Rep3/transcripts.gtf $RNA_ALIGN_DIR/HBR_Rep3.bam

stringtie -p 4 -G $RNA_REF_GTF -l UHR_Rep1 -o UHR_Rep1/transcripts.gtf $RNA_ALIGN_DIR/UHR_Rep1.bam
stringtie -p 4 -G $RNA_REF_GTF -l UHR_Rep2 -o UHR_Rep2/transcripts.gtf $RNA_ALIGN_DIR/UHR_Rep2.bam
stringtie -p 4 -G $RNA_REF_GTF -l UHR_Rep3 -o UHR_Rep3/transcripts.gtf $RNA_ALIGN_DIR/UHR_Rep3.bam
