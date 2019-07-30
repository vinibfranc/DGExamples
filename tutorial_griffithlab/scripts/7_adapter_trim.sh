echo $RNA_DATA_TRIM_DIR
mkdir -p $RNA_DATA_TRIM_DIR

echo $RNA_REFS_DIR
mkdir -p $RNA_REFS_DIR
cd $RNA_REFS_DIR
wget http://genomedata.org/rnaseq-tutorial/illumina_multiplex.fa

cd $RNA_HOME

flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz --reads2 $RNA_DATA_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz --target $RNA_DATA_TRIM_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz --reads2 $RNA_DATA_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz --target $RNA_DATA_TRIM_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz --reads2 $RNA_DATA_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz --target $RNA_DATA_TRIM_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22

flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz --reads2 $RNA_DATA_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz --target $RNA_DATA_TRIM_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz --reads2 $RNA_DATA_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz --target $RNA_DATA_TRIM_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz --reads2 $RNA_DATA_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz --target $RNA_DATA_TRIM_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22

cd $RNA_DATA_TRIM_DIR
fastqc *.fastq.gz
