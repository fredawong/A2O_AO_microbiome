# DESeq2
setwd("D:/study")

library(DESeq2)

# load data [Separated counts matrix and group infromation]
cts <- as.matrix(read.table("2.phylum.relativeAbundance.tsv", sep = "\t",header = T,
                     row.names = 1))
cts <- round(cts*10000)
coldata <- read.table("0.sampleGroup.tsv", sep = "\t",header = T,
                      row.names = 1)

# cts and coldata should be consistent in sample order.
cts <- cts[, rownames(coldata)]
(all(rownames(coldata) %in% colnames(cts))) & (all(rownames(coldata) == colnames(cts)))
# should be TRUE

# Construct DESeq data matrix
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ group)
dds

# Pre-filtering (see vignette for independent filtering)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Specifying the reference level
dds$group <- relevel(dds$group, ref = "AO")

# Differential expression analysis
dds <- DESeq(dds)
fc <- results(dds, contrast=c("group","A2O","AO"), 
                 alpha = 0.05) # FDR < 0.05
fc
summary(fc)
mcols(fc)$description

# Shrinkage of effect size (LFC estimates)
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="group_A2O_vs_AO", 
                    type="apeglm", parallel = T)
# plot volcano
xlim <- c(1e1,1e5); ylim <- c(-2,2)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="A2O vs. AO")
# plot Count difference between the most significant entry
d <- plotCounts(dds, gene=which.min(resLFC$padj), intgroup="group", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=group, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

# Output
write.csv(as.data.frame(fc), file="2.phylumFC.csv")

# Count data transformations
# Extracting transformed values (vsd or rld)
vsd <- vst(dds, blind=FALSE) # Variance stabilizing transformation
rld <- rlog(dds, blind=FALSE) # Regularized log transformation
head(assay(vsd), 3)
head(assay(rld), 3)

# Data quality assessment by sample clustering and visualization
# Heatmap of the count matrix
library("pheatmap")
# Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(rld)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
# Principal component plot of the samples
pcaData <- plotPCA(rld, intgroup="group", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw()
