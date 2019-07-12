#!/usr/bin/env Rscript

library(limma)
library(edgeR)


args <- commandArgs(TRUE)

# Get parameters for the test
countsData        = read.table(args[1],header=TRUE, row.names="gene_id", check.names=FALSE)
colData           = read.table(args[2], header=TRUE, row.names="sample", check.names=FALSE)

# Get output files
output_diff_counts   = args[3]#snakemake@output$diff_counts

#write(date(),file=output_log)
# Load counts data
#countsData = read.table(gene_counts,header=T,row.names=1,check.names=FALSE)

# Load col data with sample specifications
#colData = read.table(sample_conditions,header=T,row.names=1,check.names=FALSE)

## remove genes with 0 count
nullGenes   <- rownames(countsData[rowSums(countsData)==0,])
countsData  <- countsData[rowSums(countsData)!=0,]

## perform limma-voom
dge <- DGEList(countsData, group=colData$condition)
v   <- voom(dge, plot=FALSE)
fit <- lmFit(v)
fit <- eBayes(fit, robust=FALSE)

# writing in a file normalized counts
# FIXME add back genes with 0 counts
normalized_counts <- data.frame(id=row.names(v$E),v$E,row.names=NULL)
#write.table(normalized_counts, file=norm_counts, sep="\t", row.names=F, col.names=T, quote=F)

# Write DEGs to differentially_expressed_genes file
results <- topTable(fit, number = Inf)
results <- results[,c("AveExpr", "logFC", "P.Value","adj.P.Val")]
colnames(results) <- c("baseMean", "log2FoldChange", "pvalue", "padj")

#Add missing columns and re-order them to match DESeq2 output
results$lfcSE <- NA
results$stat <- NA
results <- results[,c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
write.table(results,file=output_diff_counts,sep="\t",quote=FALSE)
