#!/bin/bash
set -e
set -u
set -o pipefail

#echo $RNA_ALIGN_DIR
cd $RNA_ALIGN_DIR
#find *.bam -exec echo samtools index {} \; | sh

# BAM Read Counting
cd $RNA_HOME
mkdir bam_readcount
cd bam_readcount

# Create faidx indexed reference sequence file for use with mpileup
echo $RNA_REF_FASTA
samtools faidx $RNA_REF_FASTA

samtools mpileup -f $RNA_REF_FASTA -r 22:18918457-18918467 $RNA_ALIGN_DIR/UHR.bam $RNA_ALIGN_DIR/HBR.bam

# Create a bed file with some positions of interest
echo "22 38483683 38483683"
echo "22 38483683 38483683" > snvs.bed

# Count reference and variant bases at a specific position
bam-readcount -l snvs.bed -f $RNA_REF_FASTA $RNA_ALIGN_DIR/UHR.bam 2>/dev/null 1>UHR_bam-readcounts.txt
bam-readcount -l snvs.bed -f $RNA_REF_FASTA $RNA_ALIGN_DIR/HBR.bam 2>/dev/null 1>HBR_bam-readcounts.txt

# Parse the read counts for each base
cat UHR_bam-readcounts.txt | perl -ne '@data=split("\t", $_); @Adata=split(":", $data[5]); @Cdata=split(":", $data[6]); @Gdata=split(":", $data[7]); @Tdata=split(":", $data[8]); print "UHR Counts\t$data[0]\t$data[1]\tA: $Adata[1]\tC: $Cdata[1]\tT: $Tdata[1]\tG: $Gdata[1]\n";'
cat HBR_bam-readcounts.txt | perl -ne '@data=split("\t", $_); @Adata=split(":", $data[5]); @Cdata=split(":", $data[6]); @Gdata=split(":", $data[7]); @Tdata=split(":", $data[8]); print "HBR Counts\t$data[0]\t$data[1]\tA: $Adata[1]\tC: $Cdata[1]\tT: $Tdata[1]\tG: $Gdata[1]\n";'
