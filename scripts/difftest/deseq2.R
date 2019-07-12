#!/usr/bin/env Rscript

#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")

library("DESeq2")

#parallel <- FALSE
#if (snakemake@threads > 1) {
    #library("BiocParallel")
    # setup parallelization
    #register(MulticoreParam(snakemake@threads))
    #parallel <- TRUE
#}
args = commandArgs(trailingOnly=TRUE)

## Defining the code for DESEQ2 DE expression
data        = read.table(args[1],header=TRUE, row.names="gene_id", check.names=FALSE)
design      = read.table(args[2], header=TRUE, row.names="sample", check.names=FALSE)
#outdir     = paste(as.character(args[3]), "/", sep='')

output_diff_counts   = args[3]#snakemake@output$diff_counts

## remove genes with 0 count
nullGenes   <- rownames(data[rowSums(data)==0,])
data        <- data[rowSums(data)!=0,]

## Starting with DEseq count and design preparation, and performing the test
dds           = DESeqDataSetFromMatrix(countData=data,colData=design, design=~condition)

# Remove uninformative columns
dds           = dds[ rowSums(counts(dds)) > 1, ]
dds$condition = factor(dds$condition, levels=c("control","disease"))
DESeq2Table1  = dds
DE_res        = DESeq(DESeq2Table1)
#contrast     = c("condition", snakemake@params[["contrast"]])
res           = results(DE_res, parallel=parallel,pAdjustMethod = "none")

# shrink fold changes for lowly expressed genes
#res <- lfcShrink(dds,contrast=c("control","disease"), res=res)
# sort by p-value
res <- res[order(res$padj),]

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =FALSE)), by = 'row.names', sort = FALSE)

#store results
#svg(paste(outdir,"ma-plot.svg",sep=""))
#plotMA(res, ylim=c(-2,2))
#dev.off()

## FOR PCA plot
#pca_label = design[,1]

# obtain normalized counts
counts = rlog(DE_res, blind=FALSE)
#svg(paste(outdir,"PCA_results.svg",sep=""))
#plotPCA(counts, intgroup=c("condition"))
#dev.off()

write.table(resdata,file=output_diff_counts,quote=FALSE,sep="\t",row.names = FALSE)
