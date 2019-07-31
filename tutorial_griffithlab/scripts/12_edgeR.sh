
#!/bin/bash
set -e
set -u
set -o pipefail

# edgeR Analysis
cd $RNA_HOME/
mkdir -p de/htseq_counts
cd de/htseq_counts

# Create a mapping file to go from ENSG IDs (which htseq-count output) to Symbols
perl -ne 'if ($_ =~ /gene_id\s\"(ENSG\S+)\"\;/) { $id = $1; $name = undef; if ($_ =~ /gene_name\s\"(\S+)"\;/) { $name = $1; }; }; if ($id && $name) {print "$id\t$name\n";} if ($_=~/gene_id\s\"(ERCC\S+)\"/){print "$1\t$1\n";}' $RNA_REF_GTF | sort | uniq > ENSG_ID2Name.txt
head ENSG_ID2Name.txt

cut -f 1 ENSG_ID2Name.txt | sort | uniq | wc
cut -f 2 ENSG_ID2Name.txt | sort | uniq | wc
cut -f 2 ENSG_ID2Name.txt | sort | uniq -c | sort -r | head

R

