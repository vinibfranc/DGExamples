## ----style-Sweave, eval=FALSE, echo=FALSE, results="asis"----------------
#  BiocStyle::latex()

## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(
  concordance=TRUE
)

## ----Rsetup, include=FALSE------------------------------------------------------------------------
set.seed(20161000)
options(digits=2)
options(width=100)
options(browser="false")

## ----Rgetstarted----------------------------------------------------------------------------------
library(Glimma)
library(limma)
library(edgeR)
#data(lymphomaRNAseq)
data("countData.txt")
rnaseq <- lymphomaRNAseq
rnaseq$samples$group

## ----Rpreprocessing-------------------------------------------------------------------------------
rnaseq <- rnaseq[rowSums(cpm(rnaseq)>1)>=3,]
rnaseq <- calcNormFactors(rnaseq)

## ----Rmds-----------------------------------------------------------------------------------------
groups <- rnaseq$samples$group
glMDSPlot(rnaseq, groups=groups)

## ----Rvoomqw--------------------------------------------------------------------------------------
design <- model.matrix(~0+groups)
contrasts <- cbind(Smchd1null.vs.WT=c(-1,1))
vm <- voomWithQualityWeights(rnaseq, design=design)
fit <- lmFit(vm, design=design)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit)
dt <- decideTests(fit)
summary(dt)

## ----Rmd------------------------------------------------------------------------------------------
glMDPlot(fit, status=dt, counts=vm, groups=groups, side.main="Symbols")

## ----Rxy------------------------------------------------------------------------------------------
glXYPlot(x=fit$coef, y=fit$lod, xlab="logFC", ylab="logodds",
         status=dt, counts=vm, groups=groups, anno=fit$genes)

## ----Rmdsmatrix, eval=TRUE------------------------------------------------------------------------
lcpm <- cpm(rnaseq, log=TRUE, normalized.lib.sizes=TRUE)
glMDSPlot(lcpm, groups=groups)

## ----Rmdscheckinput, eval=FALSE, echo=FALSE-------------------------------------------------------
#  # Check that MDS plot runs on EList object: PASS
#  glMDSPlot(vm, groups=groups)
#  # Check that MDS plot runs on DESeqDataSet object: FAIL regarding groups
#  library(DESeq2)
#  rnaseq.deseq2 <- DESeqDataSetFromMatrix(rnaseq$counts, colData=rnaseq$samples, design=~group)
#  glMDSPlot(rnaseq.deseq2, groups=groups)

## ----Rmdsgroups-----------------------------------------------------------------------------------
groups.df <- as.data.frame(cbind(
  genotype=as.character(groups),
  lane=c(rep(4,6),3)))
groups.df
glMDSPlot(lcpm, groups=groups.df)

## ----Rmdsimple------------------------------------------------------------------------------------
glMDPlot(fit)

## ----Rmdage---------------------------------------------------------------------------------------
groups.age <- runif(ncol(rnaseq), min=3, max=12)
groups.age

## ----Rmdsamples-----------------------------------------------------------------------------------
cols <- c("yellow","blue","magenta")
sample.cols <- c("limegreen", "purple")[groups]
glMDPlot(fit, status=dt, counts=vm, groups=groups.age,
         sample.cols=sample.cols, cols=cols,
         side.ylab="logCPM", side.xlab="Age (in months)",
         side.main="Symbols", main=colnames(dt))

## ----Rmdcheckinput, eval=FALSE, echo=FALSE--------------------------------------------------------
#  # limma - PASS
#  glMDPlot(fit)
#  # edgeR - PASS
#  rnaseq.edger <- estimateDisp(rnaseq, design=design)
#  fit.edger <- exactTest(rnaseq.edger)
#  glMDPlot(fit.edger)
#  # DESeq2 - FAIL regarding annotation
#  library(DESeq2)
#  rnaseq.deseq2 <- DESeqDataSetFromMatrix(rnaseq$counts, colData=rnaseq$samples, design=~group)
#  mcols(rnaseq.deseq2) <- DataFrame(mcols(rnaseq.deseq2), rnaseq$genes)
#  fit.deseq2 <- DESeq(rnaseq.deseq2)
#  glMDPlot(fit.deseq2)

## ----Rmdanno--------------------------------------------------------------------------------------
ID <- paste(fit$genes$Symbols, fit$genes$GeneID)
DE <- c("downregulated", "notDE", "upregulated")[as.factor(dt)]
anno <- as.data.frame(cbind(ID, DE))
head(anno)
glMDPlot(fit, counts=vm, groups=groups, side.main="ID",
         anno=anno, display.columns=c("ID", "GeneName", "DE"))

## ----Rxyvolcanosimple-----------------------------------------------------------------------------
glXYPlot(x=fit$coef, y=fit$lod)

## ----Rxyvolcanomore-------------------------------------------------------------------------------
glXYPlot(x=fit$coef, y=fit$lod, xlab="logFC", ylab="logodds",
         status=dt, anno=anno, side.main="ID",
         counts=vm, groups=groups, sample.cols=sample.cols)

## ----Rarraydata-----------------------------------------------------------------------------------
data(arraydata)
arrays <- arraydata$arrays
targets <- arraydata$targets
dim(arrays)
targets

## ----Rmds.arrays----------------------------------------------------------------------------------
glMDSPlot(arrays, groups=targets[,c("Condition", "Chip", "Experiment")])

## ----Rarraysfit-----------------------------------------------------------------------------------
design <- model.matrix(~0+targets$Condition+as.factor(targets$Experiment))
contrasts <- cbind(
  DP_Ezh2KO.vs.WT=c(1,-1,0,0,0,0),
  Lum_Ezh2KO.vs.WT=c(0,0,1,-1,0,0))
fit <- lmFit(arrays, design)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit)
dt <- decideTests(fit, p.value=0.1)
summary(dt)

## ----Rarraysmd------------------------------------------------------------------------------------
sample.cols <- c("purple", "magenta", "green")[targets$Experiment]
for (COEF in 1:2) {
  glMDPlot(fit, status=dt, coef=COEF, main=colnames(fit)[COEF],
           counts=arrays, groups=targets$Condition, sample.cols=sample.cols,
           side.ylab="Log-expression", side.main="ProbeID")
}

## ----Rarraysxy------------------------------------------------------------------------------------
dt2 <- rep(0, nrow(dt))
dt2[rowSums(dt!=0)==1] <- -1
dt2[rowSums(dt!=0)==2] <- 1
table(dt2)
cols <- c("black", "grey", "red")
glXYPlot(fit$coef[,1], y=fit$coef[,2], xlab="DP", ylab="Lum",
         status=dt2, cols=cols, anno=fit$genes, side.main="ProbeID",
         counts=arrays, groups=targets$Condition, sample.cols=sample.cols,
         side.ylab="Log-expression", main="logFCs")

## ----Rmdedgeret, eval=FALSE-----------------------------------------------------------------------
#  groups <- rnaseq$samples$group
#  design <- model.matrix(~groups)
#  colnames(design) <- c("WT", "Smchd1null.vs.WT")
#  rnaseq.edger <- estimateDisp(rnaseq, design=design)
#  fit.edger <- exactTest(rnaseq.edger)
#  dt.edger <- decideTestsDGE(fit.edger)
#  glMDPlot(fit.edger, status=dt.edger, counts=rnaseq, groups=groups, transform=TRUE)

## ----Rmdedgerlrt, eval=FALSE, echo=FALSE----------------------------------------------------------
#  fit.edger <- glmFit(rnaseq.edger, design)
#  fit.edger <- glmLRT(fit.edger)
#  dt.edger <- decideTestsDGE(fit.edger)
#  glMDPlot(fit.edger, status=dt.edger, counts=rnaseq, groups=groups, transform=TRUE)

## ----Rmddeseq2, eval=FALSE------------------------------------------------------------------------
#  # BUG regarding the scale of sample expression
#  library(DESeq2)
#  rnaseq.deseq2 <- DESeqDataSetFromMatrix(
#    rnaseq$counts, colData=rnaseq$samples, design=~group)
#  mcols(rnaseq.deseq2) <- DataFrame(mcols(rnaseq.deseq2), rnaseq$genes)
#  rnaseq.deseq2 <- DESeq(rnaseq.deseq2)
#  fit.deseq2 <- results(rnaseq.deseq2, contrast=c("group", "Smchd1-null", "WT"))
#  dt.deseq2 <- as.numeric(fit.deseq2$padj<0.05)
#  glMDPlot(fit.deseq2, status=dt.deseq2, counts=rnaseq, groups=groups, transform=FALSE,
#    samples=colnames(rnaseq), anno=rnaseq$genes)

## ----Rsessinfo------------------------------------------------------------------------------------
sessionInfo()

## ----cleanup, include=FALSE-----------------------------------------------------------------------
# unlink("glimma-plots", recursive=TRUE)
