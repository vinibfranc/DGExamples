#!/bin/bash
set -e
set -u
set -o pipefail

# Use Stringtie to merge predicted transcripts from all libraries into a unified transcriptome

# Merge all 6 Stringtie results so that they will have the same set of transcripts for comparison purposes

# For reference guided mode:
cd $RNA_HOME/expression/stringtie/ref_guided/
ls -1 *Rep*/transcripts.gtf > assembly_GTF_list.txt
cat assembly_GTF_list.txt
stringtie --merge -p 8 -o stringtie_merged.gtf -G $RNA_REF_GTF assembly_GTF_list.txt

#awk '{if($3=="transcript") print}' stringtie_merged.gtf | cut -f 1,4,9 | less

# Compare reference guided transcripts to the known annotations
gffcompare -r $RNA_REF_GTF -o gffcompare stringtie_merged.gtf
cat gffcompare.stats

#awk '{if($3=="transcript") print}' gffcompare.annotated.gtf | cut -f 1,4,9 | less

# For de novo mode

cd $RNA_HOME/expression/stringtie/de_novo/
ls -1 *Rep*/transcripts.gtf > assembly_GTF_list.txt
cat assembly_GTF_list.txt
stringtie --merge -p 8 -o stringtie_merged.gtf assembly_GTF_list.txt

# Compare the de novo merged transcripts to the known annotations
gffcompare -r $RNA_REF_GTF -o gffcompare stringtie_merged.gtf
cat gffcompare.stats

