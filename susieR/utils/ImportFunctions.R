ImportVariantList <- function(VariantListPath) {
    VariantListDf <- readr::read_tsv(VariantListPath)
    VariantListDf
}

importQtlmapCovariates <- function(covariates_path){
  pc_matrix = read.table(covariates_path, check.names = F, header = T, stringsAsFactors = F)
  pc_transpose = t(pc_matrix[,-1])
  colnames(pc_transpose) = pc_matrix$SampleID
  pc_df = dplyr::mutate(as.data.frame(pc_transpose), genotype_id = rownames(pc_transpose)) %>%
    dplyr::as_tibble() %>% 
    dplyr::select(genotype_id, dplyr::everything())
  
  #Make PCA matrix
  pc_matrix = as.matrix(dplyr::select(pc_df,-genotype_id))
  rownames(pc_matrix) = pc_df$genotype_id
  return(pc_matrix)
}
importQtlmapPermutedPvalues <- function(perm_path){
  tbl = read.table(perm_path, check.names = F, header = T, stringsAsFactors = F) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(p_fdr = p.adjust(pval_beta, method = "fdr")) %>%
    dplyr::mutate(group_id = phenotype_id)
  return(tbl)
}



