#!/usr/bin/env Rscript

# force usage of conda R environment
# according to: https://github.com/conda-forge/r-rjags-feedstock/issues/6#issuecomment-504938719
.libPaths(R.home("library"))

library("data.table")
library("ssPATHS")

################## base formating functions ######################
toChar <- function(vec){
  return(as.character(unlist(vec)))
}

remove_period <- function(GeneID_vector){
 myInterestingGenes = as.character(GeneID_vector)
 myInterestingGenes=gsub("[.][0-9]*", "", myInterestingGenes)
 return(myInterestingGenes)
}

convert_df_to_se <- function(in_df, data_col_idx, annot_col_idx, sample_id_col_idx){
    
    # we must take the dataframes and converts them to a SummarizedExperiment
    out_se <- SummarizedExperiment(t(in_df[ , data_col_idx]),
                                   colData=in_df[ , annot_col_idx])
    colnames(out_se) <- in_df[,sample_id_col_idx]
    colData(out_se)$sample_id <- in_df[,sample_id_col_idx]
    
    return(out_se)
    
}

################## read in data ######################

read_patient_file <- function(curr_file){

  curr_patient = data.frame(fread(curr_file))
  patient_id = basename(curr_file)
  patient_id = stringi::stri_extract(patient_id,regex = "[^__]+")
  patient_id = stringi::stri_extract(patient_id,regex = "[^-]+")

  colnames(curr_patient) = c("gene_id", patient_id)
  curr_patient$gene_id = remove_period(curr_patient$gene_id)

  # get hypoxia genes
  hypoxia_gene_ids = get_hypoxia_genes()
  curr_patient = curr_patient[curr_patient$gene_id %in% hypoxia_gene_ids,]

  # now format it so columns are genes, rows are samples
  curr_patient_formatted = data.frame(sample_id=patient_id, t(curr_patient[,patient_id]+1))
  colnames(curr_patient_formatted) = c("sample_id", curr_patient$gene_id)
  
  # bug in ssPATHS, currently requires more than one sample, this is a work around
  if(nrow(curr_patient_formatted) == 1){
      curr_patient_formatted = rbind(curr_patient_formatted, curr_patient_formatted)
      curr_patient_formatted$sample_id[2] = "DELETE_THIS_SAMPLE"
  }
  
  # convert it to a SummarizedExperiment
  data_col_idx = grep("ENSG", colnames(curr_patient_formatted), invert=F)
  annot_col_idx = grep("ENSG", colnames(curr_patient_formatted), invert=T)
  sample_id_col_idx = grep("sample_id", colnames(curr_patient_formatted), invert=F)
  curr_patient_formatted_se = convert_df_to_se(curr_patient_formatted, data_col_idx, annot_col_idx, sample_id_col_idx)
  

  return(curr_patient_formatted_se)

}

get_tcga_training_data <- function(tcga_expr_file){

  tcga_expr_df = data.frame(fread(tcga_expr_file))

  # now only get the samples that have atleast 4 normals and SARC
  normal_freq = table(tcga_expr_df[,c("study", "is_normal")])
  normal_freq = data.frame(normal_freq)
  normal_freq = normal_freq[(normal_freq$Freq > 3 & normal_freq$is_normal==TRUE) |
                  (normal_freq$Freq > 3 & normal_freq$is_normal==FALSE),]

  normal_freq_pass = normal_freq$study[duplicated(normal_freq$study)]
  tcga_expr_df = tcga_expr_df[tcga_expr_df$study %in% c("SKCM", toChar(normal_freq_pass)),]

  # now filter for hypoxia genes

  # get hypoxia genes
  hypoxia_gene_ids = get_hypoxia_genes()
  hypoxia_gene_ids = hypoxia_gene_ids[hypoxia_gene_ids %in% colnames(tcga_expr_df)]

  # filter for gene expression cuttoff
  mean_expr = apply(tcga_expr_df[,hypoxia_gene_ids], 2, quantile, 0.25)
  mean_expr = mean_expr[mean_expr > 1000]
  genes_interest = hypoxia_gene_ids[hypoxia_gene_ids %in% names(mean_expr)]

  hypoxia_df = tcga_expr_df[,c("tcga_id", "is_normal", genes_interest)]

  # make the appropriate labels
  colnames(hypoxia_df)[1:2] = c("sample_id", "Y")
  hypoxia_df$Y = 0
  hypoxia_df$Y[tcga_expr_df$is_normal==TRUE] = 0
  hypoxia_df$Y[tcga_expr_df$is_normal==FALSE] = 1

  hypoxia_df = unique(hypoxia_df)
  
  # convert to a SummarizedExpriment
  data_col_idx = grep("ENSG", colnames(hypoxia_df), invert=F)
  annot_col_idx = grep("ENSG", colnames(hypoxia_df), invert=T)
  sample_id_col_idx = grep("sample_id", colnames(hypoxia_df), invert=F)
  hypoxia_se = convert_df_to_se(hypoxia_df, data_col_idx, annot_col_idx, sample_id_col_idx)
      

  return(hypoxia_se)

}

################## the main function ######################

get_hypoxia_score <- function(input_file, output_file, tcga_expr_file, cancer_type){

    # get patient expression
    patient_expr = read_patient_file(input_file)

    # get training data
    tcga_training_data = get_tcga_training_data(tcga_expr_file)

    # get the gene weights
    res = get_gene_weights(tcga_training_data, gene_ids=rownames(tcga_training_data), unidirectional = F)
    gene_weights = res[[1]]
    sample_scores = res[[2]]

    # calculate the score
    new_score_df_calculated = get_new_samp_score(gene_weights, patient_expr, gene_ids=rownames(patient_expr), run_normalization=TRUE)
    
    # delete the excess samples needed to go-around an ssPATHS bug
    new_score_df_calculated = subset(new_score_df_calculated, sample_id != "DELETE_THIS_SAMPLE")

    # format the output
    new_score_df_calculated$cancer_type_match = cancer_type
    new_score_df_calculated$is_normal = FALSE
    new_score_df_calculated$study = "tumor_profiler"

    write.table(new_score_df_calculated, output_file, sep="\t", quote=F, row.names=F)

}





args = commandArgs(trailingOnly=TRUE)

cancer_type = as.character(args[1])
tcga_expr_file = as.character(args[2])
patient_expr_file = as.character(args[3])
output_file = as.character(args[4])

get_hypoxia_score(patient_expr_file, output_file, tcga_expr_file, cancer_type)
