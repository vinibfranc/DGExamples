# Alignment Free Expression Estimation (Kallisto)

# Obtain transcript sequences in fasta format
cd $RNA_HOME/refs
gtf_to_fasta $RNA_REF_GTF $RNA_REF_FASTA chr22_ERCC92_transcripts.fa

cd $RNA_HOME/refs
cat chr22_ERCC92_transcripts.fa | perl -ne 'if ($_ =~/^\>\d+\s+\w+\s+(ERCC\S+)[\+\-]/){print ">$1\n"}elsif($_ =~ /\d+\s+(ENST\d+)/){print ">$1\n"}else{print $_}' > chr22_ERCC92_transcripts.clean.fa
wc -l chr22_ERCC92_transcripts*.fa

# Create a list of all transcript IDs for later use:
cd $RNA_HOME/refs
cat chr22_ERCC92_transcripts.clean.fa | grep ">" | perl -ne '$_ =~ s/\>//; print $_' | sort | uniq > transcript_id_list.txt

# Build a Kallisto transcriptome index
# Kallisto does not perform alignment or use a reference genome sequence. Instead it performs pseudoalignment to determine the compatibility of reads with targets (transcript sequences in this case). However, similar to alignment algorithms like Tophat or STAR, Kallisto requires an index to assess this compatibility efficiently and quickly.
cd $RNA_HOME/refs
mkdir -p kallisto
cd kallisto
kallisto index --index=chr22_ERCC92_transcripts_kallisto_index ../chr22_ERCC92_transcripts.clean.fa

# Generate abundance estimates for all samples using Kallisto
echo $RNA_DATA_DIR
cd $RNA_HOME/expression/
mkdir kallisto
cd kallisto

kallisto quant --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=UHR_Rep1_ERCC-Mix1 --threads=4 --plaintext $RNA_DATA_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz
kallisto quant --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=UHR_Rep2_ERCC-Mix1 --threads=4 --plaintext $RNA_DATA_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz
kallisto quant --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=UHR_Rep3_ERCC-Mix1 --threads=4 --plaintext $RNA_DATA_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz

kallisto quant --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=HBR_Rep1_ERCC-Mix2 --threads=4 --plaintext $RNA_DATA_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz
kallisto quant --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=HBR_Rep2_ERCC-Mix2 --threads=4 --plaintext $RNA_DATA_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz
kallisto quant --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=HBR_Rep3_ERCC-Mix2 --threads=4 --plaintext $RNA_DATA_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz

# Create a single TSV file that has the TPM abundance estimates for all six samples.

cd $RNA_HOME/expression/kallisto
paste */abundance.tsv | cut -f 1,2,5,10,15,20,25,30 > transcript_tpms_all_samples.tsv
ls -1 */abundance.tsv | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){print "\t$1"}' | perl -ne 'print "target_id\tlength$_\n"' > header.tsv
cat header.tsv transcript_tpms_all_samples.tsv | grep -v "tpm" > transcript_tpms_all_samples.tsv2
mv transcript_tpms_all_samples.tsv2 transcript_tpms_all_samples.tsv
rm -f header.tsv

head transcript_tpms_all_samples.tsv
tail transcript_tpms_all_samples.tsv

# Compare transcript and gene abundance estimates from Kallisto to isoform abundance estimates from StringTie and counts from HtSeq-Count
# Create the gene version of the Kallisto TPM matrix

cd $RNA_HOME/expression/kallisto
wget https://raw.githubusercontent.com/griffithlab/rnaseq_tutorial/master/scripts/kallisto_gene_matrix.pl
chmod +x kallisto_gene_matrix.pl
./kallisto_gene_matrix.pl --gtf_file=$RNA_HOME/refs/chr22_with_ERCC92.gtf  --kallisto_transcript_matrix_in=transcript_tpms_all_samples.tsv --kallisto_transcript_matrix_out=gene_tpms_all_samples.tsv

# Create a custom transcriptome database to examine a specific set of genes
cd $RNA_HOME/refs
grep rRNA $RNA_REF_GTF > genes_chr22_ERCC92_rRNA.gtf 
gtf_to_fasta genes_chr22_ERCC92_rRNA.gtf $RNA_REF_FASTA chr22_rRNA_transcripts.fa
cat chr22_rRNA_transcripts.fa | perl -ne 'if ($_ =~/^\>\d+\s+\w+\s+(ERCC\S+)[\+\-]/){print ">$1\n"}elsif($_ =~ /\d+\s+(ENST\d+)/){print ">$1\n"}else{print $_}' > chr22_rRNA_transcripts.clean.fa
cat chr22_rRNA_transcripts.clean.fa

cd $RNA_HOME/refs/kallisto
kallisto index --index=chr22_rRNA_transcripts_kallisto_index ../chr22_rRNA_transcripts.clean.fa
