## This analysis was based on: http://www.nathalievialaneix.eu/doc/html/solution_edgeR-tomato-withcode.html

## Install biconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

## Install packages needed
BiocManager::install("edgeR")
BiocManager::install("limma")
BiocManager::install("mixOmics")
BiocManager::install("DESeq2")
BiocManager::install("HTSFilter")

## Loading required package
library(edgeR)
library(limma)
library(RColorBrewer)
library(mixOmics)
library(HTSFilter)

## Get working dir 
directory <- getwd()
dir(directory)

## Load read counts per sample file and the experiment design
rawCountTable <- read.table("countData.txt", header=TRUE, sep="\t", row.names=1)
sampleInfo <- read.table("design.csv", header=TRUE, sep=",", row.names=1)

## Look the count data 
head(rawCountTable)

## Count genes included in count file
nrow(rawCountTable)

## Create a DGEList data object containing counts and sample info
dgeFull <- DGEList(rawCountTable, group=sampleInfo$condition)
dgeFull

## Extract pseudo-counts (ie log2(K+1)) for data exploration and quality assessment
pseudoCounts <- log2(dgeFull$counts+1)
head(pseudoCounts)

## Plot a histogram for pseudo-counts in Cond.WT.Rep.1
hist(pseudoCounts[,"Cond.WT.Rep.1"])

## Plot a boxplot for pseudo-counts
boxplot(pseudoCounts, col="gray", las=3)

## MA-plot between WT or Mt samples

par(mfrow=c(1,2))
## WT1 vs WT2
# A values
avalues <- (pseudoCounts[,1] + pseudoCounts[,2])/2
# M values
mvalues <- (pseudoCounts[,1] - pseudoCounts[,2])
plot(avalues, mvalues, xlab="A", ylab="M", pch=19, main="treated")
abline(h=0, col="red")
## Mt1 vs Mt2
# A values
avalues <- (pseudoCounts[,4] + pseudoCounts[,5])/2
# M values
mvalues <- (pseudoCounts[,4] - pseudoCounts[,5])
plot(avalues, mvalues, xlab="A", ylab="M", pch=19, main="control")
abline(h=0, col="red")

## MDS for pseudo-counts (using limma package)
plotMDS(pseudoCounts)

## Heatmap for pseudo-counts (using mixOmics package)
sampleDists <- as.matrix(dist(t(pseudoCounts)))
sampleDists

cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(16)
cim(sampleDists, color=cimColor, symkey=FALSE)

## Remove genes with zero counts for all samples
dgeFull <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ],
                   group=dgeFull$samples$group)
head(dgeFull$counts)

# Estimate the normalization factors
dgeFull <- calcNormFactors(dgeFull, method="TMM")
dgeFull$samples

head(dgeFull$counts)

## Find the normalized counts and use the log2-transformation to inspect them with boxplots and a MDS
eff.lib.size <- dgeFull$samples$lib.size*dgeFull$samples$norm.factors
normCounts <- cpm(dgeFull)
pseudoNormCounts <- log2(normCounts + 1)
boxplot(pseudoNormCounts, col="gray", las=3)

plotMDS(pseudoNormCounts)

## Estimate common and tagwise dispersion
dgeFull <- estimateCommonDisp(dgeFull)
dgeFull <- estimateTagwiseDisp(dgeFull)
dgeFull

## Perform an exact test for the difference in expression between the two conditions “WT” and “Mt”
dgeTest <- exactTest(dgeFull)
dgeTest

## Independent filtering to remove low expressed genes
filtData <- HTSFilter(dgeFull)$filteredData

dgeTestFilt <- exactTest(filtData)
dgeTestFilt

## Diagnostic plot for multiple testing

## Plot an histogram of unadjusted p-values
hist(dgeTest$table[,"PValue"], breaks=50)

## Plot an histogram of unadjusted p-values after filtering
hist(dgeTestFilt$table[,"PValue"], breaks=50)

## Inspecting the results

## Extract a summary of the differential expression statistics
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table))
head(resNoFilt)

resFilt <- topTags(dgeTestFilt, n=nrow(dgeTest$table))
head(resFilt)

## Compare the number of differentially expressed genes with and without filtering (risk: 1%)

## Before independent filtering
sum(resNoFilt$table$FDR < 0.01)

## After independent filtering
sum(resFilt$table$FDR < 0.01)

## Extract and sort differentially expressed genes

## Downregulated genes
sigDownReg <- resFilt$table[resFilt$table$FDR<0.01,]
sigDownReg <- sigDownReg[order(sigDownReg$logFC),]
head(sigDownReg)

## Upregulated genes
sigUpReg <- sigDownReg[order(sigDownReg$logFC, decreasing=TRUE),]
head(sigUpReg)

## Write the results in csv files
write.csv(sigDownReg, file="sigDownReg_tomato.csv")
write.csv(sigUpReg, file="sigUpReg_tomato.csv")

## Interpreting the DE analysis results

## Create a MA plot with 1% differentially expressed genes
plotSmear(dgeTestFilt,
          de.tags = rownames(resFilt$table)[which(resFilt$table$FDR<0.01)])

## Create a Volcano plot
volcanoData <- cbind(resFilt$table$logFC, -log10(resFilt$table$FDR))
colnames(volcanoData) <- c("logFC", "negLogPval")
head(volcanoData)

plot(volcanoData, pch=19)

## Transform the normalized counts in log-counts-per-million
y <- cpm(dgeFull, log=TRUE, prior.count = 1)
head(y)

## Select 1% differentially expressed genes and produce a heatmap
selY <- y[rownames(resFilt$table)[resFilt$table$FDR<0.01 & 
                                    abs(resFilt$table$logFC)>1.5],]
head(selY)

cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)[255:1]
finalHM <- cim(t(selY), color=cimColor, symkey=FALSE)

## Gene clustering details
plot(finalHM$ddc, leaflab="none")
abline(h=10, lwd=2, col="pink")

## Using this dendrogram, we might want to cut the tree at level h=10
geneClust <- cutree(as.hclust(finalHM$ddc), h=10)
head(geneClust)

## Get the nnumber of clusters
length(unique(geneClust))

## Get the genes in cluster 1
names(which(geneClust==1))
