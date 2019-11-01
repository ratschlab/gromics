#!/usr/bin/env python
import sys

"""
Usage: python genename_to_HUGO.py <resources/tcga_boxplot/genes_of_interest.txt> <{sample}.genes.results.normalized.header.txt> <{sample}.genes.results.normalized.genes_filtered.txt> 

Purpose: Changes the ENSEMBL geneID to HUGO symbol.
    Usage scenario: This script is used in RNASeq snakemake pipeline. The 'rule filter_genes' in rules/tcga_boxplot.smk uses this script
    Usage example:   genename_to_HUGO.py <genes_of_interest.txt> <OB225.genes.results.normalized.header.txt>  > OB225.genes.results.normalized.genes_filtered.txt
    ---------------------       -----------------------------------------       -------------------------------------------------
    genes_of_interest.txt   |   OB225.genes.results.normalized.header.txt   |   OB225.genes.results.normalized.genes_filtered.txt
    ---------------------       -----------------------------------------       -------------------------------------------------
    AKT1,ENSG00000142208        gene_id OB225                                   gene_id OB225
                                ENSG00000142208.11      1067.7043               AKT1    1067.7043
Author: Tinu Thomas
Date: July, 2019
"""


def main():
    if len(sys.argv) < 3 :
        sys.stderr.write('Usage: %s <genes_of_interest.txt> <{sample}.genes.results.normalized.header.txt> <{sample}.genes.results.normalized.genes_filtered.txt>\n' % sys.argv[0])
        sys.exit(1)
    genes={}

    # Reading genes of interest file and storing genes in the dictionary 'genes'
    with open(sys.argv[1],'r') as fgene:
        for line in fgene:
            line=line.strip().split(',')
            genes[line[1]]=line[0] # stored HUGO symbol 

    # Reading the normalized data file and extracting information for genes of interest
    with open(sys.argv[2],'r') as fnorm:
        with open(sys.argv[3], 'w') as fout:
            for data in fnorm:
                data=data.strip().split('\t')
                if data[0] == 'gene_id':
                    fout.write('\t'.join(data))
                    fout.write('\n')
                else:
                    if data[0].split('.')[0] in genes:
                        fout.write('\t'.join([ genes[data[0].split('.')[0]] ] + data[1:] )) # changes the ENSEMBL geneID in the normalized file to HUGO symbol
                        fout.write('\n')
    return 0

if __name__ == "__main__":
    sys.exit(main())
