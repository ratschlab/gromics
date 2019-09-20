#!/usr/bin/env Rscript

# force usage of conda R environment
.libPaths(R.home("library"))
## Command line script for generating individual expression context figures
##
##
## Manual
## ------
##
## In command line type:
## Rscript --vanilla parse_rna_tcga_pat_data.r <TCGA_normalized> <sample_normalized> <out>
##
##  <TCGA_normalized>   - the TCGA reference cohort normalized count data
##  <sample_normalized> - the RSEM normalized sample counts
##  <out> 		- name of the out file
###############################################################################


## Retrieve command line arguments (all character)
args <- commandArgs(trailingOnly = TRUE)

## Read command line input.
TCGA <- read.table(args[1], stringsAsFactor=F, header=T, row.names=1)
sample <- read.table(args[2],stringsAsFactor=F, header=T)
out <- args[3]

## Initialize the result table.
data <-c("gene", "lower_whisker","lower_quartile", "median","upper_quartile","higher_whisker","patient_sample","normalized_expression","normalized_expression_text", "percentile")
## Start looping through the patients data.
print(length(sample[,1]))
for (i in 1:length(sample[,1])){

	## Read in the current gene, the patients value and the cohort values.
    gene<-sample[i,1]
	name<-(strsplit(gene, split="\\|")[[1]][1])
	#name<-(gene)
	pat<- c(sample[i,2])
	cohort <- t(TCGA[name,])
	## Calculate the percentile corresponding to the patient values.
	quant <- quantile(cohort,probs = seq(0, 1, by= 0.01), na.rm=TRUE)
	percentile<- max(1, which(pat >=quant))-1
    #if (pat==0 & all(quant==0)){percentile <- "NA"}
	#if (pat==0 & quant[1]>0){percentile <- "1"}
	#if (pat>0 & all(quant==0)){percentile <- "100"}	
	if (is.null(pat) & all(is.null(quant))){percentile <- "NA"}
	if (is.null(pat) & quant[1]>0){percentile <- "1"}
	if (!is.null(pat) & all(is.null(quant))){percentile <- "100"}	
	## Calculate in which part of a boxplot the patient value would be located. 
	## 	0=Below lower whisker, 1=Within lower whisker, 2=Within lower part of box.
	## 	3=Within upper part of box, 4=Within upper whisker, 5=Above upper whisker.
	stats <- boxplot(cohort, horizontal = T, plot =F)$stats
	pat.location <- max(0, which(pat >stats))
	if (pat.location == 0){
	text <- "under_expression"} 
	else if (pat.location == 1){
	text <- "low_expression"}
	else if (pat.location == 4){
        text <- "high_expression"}
	else if (pat.location == 5){
        text <- "over_expression"}
	else if (pat.location == 2 | pat.location == 3){
	text <- "normal"}
	else {
	print("Patient value does not fall in any expected range")}	
	## Append results to rest of data.
	data<-rbind(data, c(as.vector(name),as.vector(stats),pat, as.vector(pat.location),as.vector(text),as.vector(percentile)))
	}

## Define the first row as new header & then delete the first row.
colnames(data) = data[1, ] 
data = data[-1, ] 

## Write the results.
write.table(data, file=out,sep="\t",quote=FALSE,row.names=FALSE)

