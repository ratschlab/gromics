#!/usr/bin/env Rscript

# force usage of conda R environment
.libPaths(R.home("library"))

## Command line script for generating individual expression context figures
##
## Author: Nora Toussaint (based on Daniel Stekhoven's code), adjusted by Anne Richter for gene list instead of single gene (Mai 2018)
##
## Manual
## ------
##
## In command line type:
## Rscript --vanilla variant_expression_context.TCGA.R <outdir> <TCGA_normalized> <sample_normalized> <gene_list> <sample_name>
##`
##  <outdir>            - path to output directory
##  <TCGA_normalized>   - the TCGA reference cohort normalized count data
##  <sample_normalized> - the RSEM normalized sample counts
##  <gene_list>         - table with gene names to be plotted in first column
##  <sample_name>       - name of sample that is being analysed
##
## A PNG <outdir/<sample_name>_geneName_boxplot_tcga_cohort.png> will be generated for each gene in the list; with the boxplot in it.
###############################################################################

# # get example data
# require(readr)
# dat <- read_delim("~/gitlab/util/plot/Tumor_P2_RNA.htseq_counts.txt_TPM_HGNC.tsv.z_score_table.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
# pat <- dat$TPM[which(dat$hgnc_symbol=="DDR2")]

## Retrieve command line arguments (all character)
args <- commandArgs(trailingOnly = TRUE)
gene_table <- read.table(args[4], sep = "\t", header = T)
gene_list <- gene_table[[1]]
samplename <- args[5]

###   gene_symbol <- args[4]

## load reference data
TCGA <- read.table(args[2], header=T, row.names=1)
sample <- read.table(args[3], header=T, row.names=1)

for (gene_symbol in gene_list){
	## get full gene name
	gene <- rownames(TCGA)[grep(paste0("^", gene_symbol, "\\|"), rownames(TCGA))]
    name<-(strsplit(gene, split="\\|")[[1]][1])
    #print(gene)
    #print(name)
    #print("----")
	
    ## get relevant data
	cohort = t(TCGA[gene,])
	pat = sample[name,1]

	## generate boxplot to get stats
	bp <- boxplot(cohort, horizontal = T, plot = F)
	stats <- bp$stats

	## assign shape and color to patient expression level of gene
	## green/round 25-75%, orange/square inside whiskers, red/triangle outside whiskers
	#pat <- sample[gene,1]
	#pat.location <- max(0, which(pat >=stats))
	pat.location <- max(0, which(pat>=stats[-3]))+1
	tab.key <- data.frame(
	  location=1:5, 
	  pch=c(25, 23, 22, 23, 24), # use example(points) for other shapes
	  col=colorRampPalette(c("blue", "white", "red"))(5),
	  stringsAsFactors = F)
	pat.shape <- tab.key$pch[pat.location]
	pat.col <- tab.key$col[pat.location]


	## count number of outliers
	no.outliers <- length(bp$out)
	prop.outliers <- round(100*no.outliers/length(cohort))

	## plot final figure
	plotname <- paste(args[1], "/", args[5], "_", gene_symbol, "_boxplot_tcga_cohort.png", sep="")
    png(plotname, width = 600, height = 120)
    # set margins to minimum and remove box
	par(mar=c(0, 0, 0, 0), mai=c(0, 0, 0, 0), bty="n")
	# plot boxplot, no axis, no outliers, zoom in on box
	boxplot(cohort, horizontal = T, ylim=c(0.9*min(pat, stats[1]), 1.1*max(pat, stats[5])), outline = F, boxfill="grey90", xaxt="n", lwd = 8, boxwex = 1)
	# add patient expression level for given gene
	points(pat, 1, pch = pat.shape, col = "black", cex = 8, bg = pat.col, lwd=4)
	# (optional) print the number/proportion of outliers
	# if (print.out)
	#   text(stats[5], 0.7, paste0(no.outliers, " (", prop.outliers, "%) ", intToUtf8(9658)), adj = 0.82, cex = 1.5)
	dev.off()

	cat(no.outliers, "(", prop.outliers, "%) samples with outside whisker expression level\n")
}

