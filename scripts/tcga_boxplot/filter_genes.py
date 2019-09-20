#!/usr/bin/env python
import sys

"""
Usage: python filter_genes.py <resources/tcga_boxplot/genes_of_interest.txt> <{sample}.genes.results.normalized.header.txt>  > {sample}.genes.results.normalized.genes_filtered.txt 
Purpose: Changes the ENSEMBL geneID in the file {sample}.genes.results.normalized.header.txt to HUGO symbol.
The resources/tcga_boxplot/genes_of_interest.txt has the 'HUGO_symbol->ENSEMBL_ID' for our genes of interest
The 'rule filter_genes' in rules/tcga_boxplot.smk uses this script
Author: Tinu Thomas
Date: July, 2019
"""
genes={}

# Reading genes of interest file and storing genes in the dictionary 'genes'
with open(sys.argv[1],'r') as fgene:
    for line in fgene:
        line=line.strip().split(',')
        genes[line[1]]=line[0] # stored HUGO symbol 

# Reading the normalized data file and extracting information for genes of interest
with open(sys.argv[2],'r') as fnorm:
    for data in fnorm:
        data=data.strip().split('\t')
        if data[0] == 'gene_id':
            print('\t'.join(data))
        else:
            if data[0].split('.')[0] in genes:
                print('\t'.join([ genes[data[0].split('.')[0]] ] + data[1:] )) # changes the ENSEMBL geneID in the normalized file to HUGO symbol
