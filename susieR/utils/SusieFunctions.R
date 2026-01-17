filterMAF <- function(genotype_matrix,MAF_threshold = 0,variant_list = NULL) {
    
    if (is.null(variant_list)) {
    # computes MAF across genotype matrix
    # and filters to those above the MAF threshold
    MAF_calculations <- genotype_matrix %>% 
        t() %>% 
        data.frame() %>% 
        summarize(across(everything(),~sum(.)/(dplyr::n()*2))) %>% 
        t() %>% 
        data.frame() %>% 
        dplyr::rename('AF' = 1) %>% 
        mutate(MAF = case_when(AF > .5 ~ 1 - AF,TRUE ~ AF))  %>%
        filter(MAF > MAF_threshold)
    filtered_genotype_matrix <- genotype_matrix[rownames(MAF_calculations),]
    } else {
        # uses variant list to filter genotype matrix 
        variant_df <- ImportVariantList(variant_list)
        valid_variants <- intersect(rownames(genotype_matrix),variant_list)
        filtered_genotype_matrix <- genotype_matrix[valid_variants,]
    }
    return(filtered_genotype_matrix) 
}



finemapPhenotype <- function(phenotype_id, se, genotype_file, covariates, cis_distance,MAF = 0,variant_list = NULL){
  message("Processing phenotype: ", phenotype_id)
  
  #Extract phenotype from SE
  gene_vector = eQTLUtils::extractPhentypeFromSE(phenotype_id, se, "counts") %>%
    dplyr::mutate(phenotype_value_std = qnorm((rank(phenotype_value, na.last = "keep") - 0.5) / sum(!is.na(phenotype_value))))
  selected_phenotype = phenotype_id
  gene_meta = dplyr::filter(SummarizedExperiment::rowData(se) %>% as.data.frame(), phenotype_id == selected_phenotype)

  #Rearrange samples in the covariates matrix
  covariates_matrix = cbind(covariates[gene_vector$genotype_id,], 1)
  
  #Import genotype matrix
  genotype_matrix = eQTLUtils::extractGenotypeMatrixFromDosage(
    chr = gene_meta$chromosome, 
    start = gene_meta$phenotype_pos - cis_distance, 
    end = gene_meta$phenotype_pos + cis_distance, 
    dosage_file = genotype_file) %>% 
    filterMAF(MAF_threshold = MAF,variant_list = variant_list) 

  #Residualise gene expression and genotype matrix
  hat = diag(nrow(covariates_matrix)) - covariates_matrix %*% solve(crossprod(covariates_matrix)) %*% t(covariates_matrix)
  expression_vector = hat %*% gene_vector$phenotype_value_std
  names(expression_vector) = gene_vector$genotype_id
  
  # subset genotype matrix to individuals that are in 
  # expression vector
  gt_matrix = genotype_matrix[,names(expression_vector)]
  
  #Exclude variants with no alternative alleles
  gt_matrix = gt_matrix[rowSums(round(gt_matrix,0), na.rm = TRUE) != 0,]
  
  #Replace missing values with row means
  gt_matrix = t(gt_matrix) %>% zoo::na.aggregate() %>% t()

  #Standardise genotypes
  gt_std = t(gt_matrix - apply(gt_matrix, 1, mean))
  gt_hat = hat %*% gt_std
  
  # Fit finemapping model
  fitted <- susieR::susie(gt_hat, expression_vector,
                          L = 10,
                          estimate_residual_variance = TRUE, 
                          estimate_prior_variance = TRUE,
                          scaled_prior_variance = 0.1,
                          verbose = TRUE,
                          compute_univariate_zscore = TRUE,
                          min_abs_corr = 0.5)
  fitted$variant_id = rownames(gt_matrix)
  return(fitted)
}


extractResults <- function(susie_object){
  credible_sets = susie_object$sets$cs
  cs_list = list()
  susie_object$sets$purity = dplyr::as_tibble(susie_object$sets$purity) %>%
    dplyr::mutate(
      cs_id = rownames(susie_object$sets$purity),
      cs_size = NA,
      cs_log10bf = NA,
      overlapped = NA
    )
  added_variants = c()
  for (index in seq_along(credible_sets)){
    cs_variants = credible_sets[[index]]
    cs_id = susie_object$sets$cs_index[[index]]

    is_overlapped = any(cs_variants %in% added_variants)
    susie_object$sets$purity$overlapped[index] = is_overlapped
    susie_object$sets$purity$cs_size[index] = length(cs_variants)
    susie_object$sets$purity$cs_log10bf[index] = log10(exp(susie_object$lbf[cs_id]))
    if (!is_overlapped) {
      cs_list[[index]] = dplyr::tibble(cs_id = paste0("L", cs_id),
                                       variant_id = susie_object$variant_id[cs_variants])
      added_variants = append(added_variants, cs_variants)
    }
  }
  df = purrr::map_df(cs_list, identity)

  #Extract purity values for all sets
  purity_res = susie_object$sets$purity

  #Sometimes all the PIP values are 0 and there are no purity values, then skip this step
  if(nrow(purity_res) > 0){
    purity_df = dplyr::as_tibble(purity_res) %>%
      dplyr::filter(!overlapped) %>%
      dplyr::mutate(
        cs_avg_r2 = mean.abs.corr^2,
        cs_min_r2 = min.abs.corr^2,
        low_purity = min.abs.corr < 0.5
      )  %>%
      dplyr::select(cs_id, cs_log10bf, cs_avg_r2, cs_min_r2, cs_size, low_purity) 
  } else{
    purity_df = dplyr::tibble()
  }
  
  #Extract betas and standard errors and lbf_variables
  mean_vec = susieR::susie_get_posterior_mean(susie_object)
  sd_vec = susieR::susie_get_posterior_sd(susie_object)
  
  #Extract matrices
  alpha_mat = t(susie_object$alpha)
  colnames(alpha_mat) = paste0("alpha", seq(ncol(alpha_mat)))
  
  mu_mat = t(susie_object$mu)
  colnames(mu_mat) = paste0("mu_", seq(ncol(mu_mat)))
  
  mu2_mat = t(susie_object$mu2)
  colnames(mu2_mat) = paste0("mu2_", seq(ncol(mu2_mat)))
  
  lbf_variable_mat = t(susie_object$lbf_variable)
  colnames(lbf_variable_mat) = paste0("lbf_variable", seq(ncol(lbf_variable_mat)))
  posterior_df = dplyr::tibble(variant_id = rownames(alpha_mat), 
                               #pip = susie_object$pip,
                               z = susie_object$z[,1],
                               posterior_mean = mean_vec, 
                               posterior_sd = sd_vec,
                               X_column_scale_factors = susie_object$X_column_scale_factors) %>%
                 dplyr::bind_cols(purrr::map(list(alpha_mat, mu_mat, mu2_mat), dplyr::as_tibble))
  lbf_df = dplyr::tibble(variant_id = rownames(lbf_variable_mat)) %>%
    dplyr::bind_cols(dplyr::as_tibble(lbf_variable_mat))

  if(nrow(df) > 0 & nrow(purity_df) > 0 & ncol(lbf_df) > 10){ #ncol(lbf_df) <= 10 only if the number of variants in the region is < 10
    cs_df = purity_df
    variant_df = dplyr::left_join(posterior_df, df, by = "variant_id") %>%
      dplyr::left_join(cs_df, by = "cs_id")
  } else{
    cs_df = NULL
    variant_df = NULL
    lbf_df = NULL
  }

  return(list(cs_df = cs_df, variant_df = variant_df, lbf_df = lbf_df))
}

extractPipsFromVariantDf <- function(variant_df){
  alpha_df = dplyr::select(variant_df, phenotype_id, variant_id, cs_id, cs_index, alpha1:alpha10)
  #Rename alpha1:alpha10 to L1:L10
  colnames(alpha_df) = c("phenotype_id", "variant_id", "cs_id","cs_index", "L1","L2","L3","L4", "L5", "L6", "L7", "L8", "L9", "L10")
  
  pip_df = dplyr::filter(alpha_df, !is.na(cs_index)) %>% 
    tidyr::pivot_longer(L1:L10, values_to = "pip") %>% 
    dplyr::filter(cs_index == name) %>% 
    dplyr::select(-name)
  
  return(pip_df)
}

make_connected_components_from_cs <- function(susie_all_df, z_threshold = 3, cs_size_threshold = 10) {
  # Filter the credible sets by Z-score and size
  susie_filt_all <- susie_all_df %>%
    dplyr::group_by(group_id) %>%
    dplyr::mutate(max_abs_z = max(abs(z))) %>%
    dplyr::filter(max_abs_z > z_threshold, cs_size < cs_size_threshold) %>%
    dplyr::ungroup()
  
  susie_highest_pip_per_cc <- data.frame()
  
  uniq_groups = susie_filt_all$group_id %>% base::unique()
  # make the ranges object in order to find overlaps
  for (uniq_group in uniq_groups) {
    susie_filt <- susie_filt_all %>% dplyr::filter(group_id == uniq_group)
    message("Processing CC of group_id: ", uniq_group)

    cs_ranges = GenomicRanges::GRanges(
      seqnames = susie_filt$chromosome,
      ranges = IRanges::IRanges(start = susie_filt$position, end = susie_filt$position),
      strand = "*",
      mcols = data.frame(molecular_trait_id = susie_filt$molecular_trait_id, variant_id = susie_filt$variant, group_id = susie_filt$group_id)
    )
    
    # find overlaps and remove the duplicated
    olaps <- GenomicRanges::findOverlaps(cs_ranges, cs_ranges, ignore.strand = TRUE) %>%
      GenomicRanges::as.data.frame() %>%
      dplyr::filter(queryHits <= subjectHits)
    
    # change variant sharing into credible set sharing 
    # not to have multiple connected components of variants but credible sets
    olaps <- olaps %>% dplyr::mutate(cs_mol_1 = cs_ranges$mcols.molecular_trait_id[queryHits], cs_mol_2 = cs_ranges$mcols.molecular_trait_id[subjectHits])
    edge_list <- olaps %>% dplyr::select(cs_mol_1, cs_mol_2) %>% BiocGenerics::unique() %>% base::as.matrix()
    
    # make the graph of connected components
    g <- igraph::graph_from_edgelist(edge_list, directed = F)
    g_cc <- igraph::components(g)
    
    # turn connected components graph into data frame
    cc_df <- data.frame(cc_membership_no = g_cc$membership, 
                        molecular_trait_id = g_cc$membership %>% names()) 
    
    susie_highest_pip_per_cc_temp <- susie_filt %>% 
      dplyr::left_join(cc_df, by = "molecular_trait_id") %>% 
      dplyr::group_by(cc_membership_no) %>% 
      dplyr::arrange(-pip) %>% 
      dplyr::slice(1) %>% 
      dplyr::ungroup()
    
    susie_highest_pip_per_cc <- susie_highest_pip_per_cc %>% base::rbind(susie_highest_pip_per_cc_temp)
  }
  
  return(susie_highest_pip_per_cc)
}

splitIntoBatches <- function(n, batch_size){
  n_batches = ceiling(n/batch_size)
  batch_ids = rep(seq(1:n_batches), each = batch_size)[1:n]
  return(batch_ids)
}

splitIntoChunks <- function(chunk_number, n_chunks, n_total){
  chunk_size = max(1,floor(n_total/(n_chunks)))
  batches = splitIntoBatches(n_total,chunk_size)
  batches[batches > n_chunks] = seq(from = 1, to = length(batches[batches > n_chunks]))
  selected_batch = batches == chunk_number
  return(selected_batch)
}
