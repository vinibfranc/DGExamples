#!/bin/bash
set -e
set -u
set -o pipefail

cd $RNA_HOME/data
fastqc *.fastq.gz

cd $RNA_HOME/data
multiqc .
