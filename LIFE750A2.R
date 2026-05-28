source("install_bioc_packages.R")

library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(gplots)

# Read in counts and metadata files, sep="\t" for TSV
counts_data <- read.table("gene_counts.tsv", header=TRUE, row.names=1, sep="\t")
head(counts_data)
colData <- read.table("sample_metadata.tsv", header=TRUE, row.names=1, sep="\t")

# Column names in counts match the row names in metadata
all(colnames(counts_data) %in% rownames(colData))
all(colnames(counts_data) == rownames(colData))
# Manually match the column name with row names in metadata
counts_data <- counts_data[, rownames(colData)]
# Check again
all(colnames(counts_data) == rownames(colData))

# Construct the DESeqdataset object
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ condition)
dds

# Pre-filtering: remove genes less than 10 reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

# Set the factor level (The Baseline)
dds$Condition <- relevel(dds$condition, ref = "normal")

# Run the DESeq analysis
dds <- DESeq(dds)

# Extract results and order by most significant p-value
res <- results(dds)
res <- res[order(res$padj), ]

summary(res)

# Stricter 0.01 FDR threshold
res0.01 <- results(dds, alpha = 0.01)

res0.01 <- res0.01[order(res0.01$padj), ]
summary(res0.01)

# contrasts, [Metadata_Column]_[Test_Group]_vs_[Baseline_Group]
resultsNames(dds)

# Force normal to be the baseline
dds$condition <- relevel(dds$condition, ref = "normal")

# Re run
dds <- DESeq(dds)

# Extract results and order by most significant p-value
res <- results(dds)
res <- res[order(res$padj), ]

summary(res)

# Stricter 0.01 FDR threshold
res0.01 <- results(dds, alpha = 0.01)

res0.01 <- res0.01[order(res0.01$padj), ]
summary(res0.01)

# Check the name again
resultsNames(dds)

# MA Plot to visualize up vs downregulation
plotMA(res0.01)

# Merge with the normalized count data 
resdata <- merge(as.data.frame(res0.01), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene" 
head(resdata)

# Write the results to a text file
write.table(resdata, file="diffexpr_results.txt", sep="\t", quote=F)

# Filter for genes where the adjusted p-value < 0.01 and have no NA value
sig_genes <- resdata[!is.na(resdata$padj) & resdata$padj < 0.01, ]
sig_genes

# Look directly at the 0.01 strict DESeq2 object as only able to extract gene ID from here 
sig_res <- res0.01[!is.na(res0.01$padj) & res0.01$padj < 0.01, ]
sig_res

# Find the most down and upregulated gene
most_down_gene <- rownames(sig_res)[which.min(sig_res$log2FoldChange)]
most_up_gene <- rownames(sig_res)[which.max(sig_res$log2FoldChange)]
most_down_gene
most_up_gene

# plotCounts 
plotCounts(dds, gene=most_down_gene, intgroup="condition", main="Most Downregulated Gene")
plotCounts(dds, gene=most_up_gene, intgroup="condition", main="Most Upregulated Gene")

# Single most significantly changed gene 
plotCounts(dds, gene=rownames(res0.01)[1], intgroup="condition") 

# Info about variable and test
mcols(res0.01)$description

# Plot of adjusted p-value
hist(res$padj, breaks=50, col="grey", main="Histogram of Adjusted P-values")

# Perform regularized log transformation for heatmaps/PCA 
rld <- rlogTransformation(dds) 
head(assay(rld))
hist(assay(rld))

# Sample Distance Heatmap
sampleDists <- as.matrix(dist(t(assay(rld)))) 
mycols <- brewer.pal(8,"Dark2")[1:length(unique(colData$condition))]
heatmap.2(as.matrix(sampleDists), key=F, trace="none", col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[colData$condition], RowSideColors=mycols[colData$condition],
          margin=c(10, 10), main="Sample Distance Matrix") 

# Check similarity or differenecs by PCA plot 
plotPCA(rld, intgroup=c("Condition")) 
