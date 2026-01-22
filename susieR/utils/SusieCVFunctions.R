RunFoldCV <- function(Metadata,
                      CovariatesMatrix,
                      CVObj,
                      ExpressionMatrix,
                      GenotypeMatrix,
                      GeneID,
                      PhenotypeMeta,
                      CisDistance,
                      RegionDf,
                      Fold,
                      VariantList = NULL
                      ) {
    FoldMetadata <-  Metadata %>% mutate(qtl_group = case_when(fold == k ~ 'Test',TRUE ~ 'Train'))
    TestSamples <- FoldMetadata %>% filter(qtl_group == 'Test') %>% mutate(sample_id = as.character(sample_id)) 
    TrainSamples <- FoldMetadata %>% filter(qtl_group == 'Train') %>% mutate(sample_id = as.character(sample_id))
    TrainCovariateSet <- GetFoldTrainPCs(CVObj,k) %>% MergeMolecularGeneticPCs(CovariatesMatrix)
    TestCovariateSet <- GetFoldTestPCs(CVObj,k) %>% MergeMolecularGeneticPCs(CovariatesMatrix)
    TrainBedDf <- expression_matrix %>% select(1,2,3,4,all_of(TrainSamples$sample_id))
    TestBedDf <- expression_matrix %>% select(1,2,3,4,all_of(TestSamples$sample_id))

    TrainSe <- eQTLUtils::makeSummarizedExperimentFromCountMatrix(assay = TrainBedDf %>% 
                                                                        select(-1,-2,-3) , 
                                                                    row_data = PhenotypeMeta, 
                                                                    col_data = TrainSamples, 
                                                                    quant_method = "gene_counts",
                                                                    reformat = FALSE)
    TestSe <- eQTLUtils::makeSummarizedExperimentFromCountMatrix(assay = TestBedDf %>% 
                                                                        select(-1,-2,-3) , 
                                                                    row_data = PhenotypeMeta, 
                                                                    col_data = TestSamples, 
                                                                    quant_method = "gene_counts",
                                                                    reformat = FALSE)

    TrainGeneVector <- eQTLUtils::extractPhentypeFromSE(GeneID, TrainSe, "counts")
    TestGeneVector <- eQTLUtils::extractPhentypeFromSE(GeneID, TestSe, "counts")
    TrainCovarModel <- EstimateBetaHat(TrainGeneVector$phenotype_value,TrainCovariateSet)     
    
    selected_phenotype <- GeneID 
    selected_qtl_group <- eQTLUtils::subsetSEByColumnValue(TrainSe, "qtl_group",'Train')
    selected_phenotypes = phenotype_list %>%
      dplyr::filter(group_id %in% selected_group_ids) %>%
      dplyr::pull(phenotype_id) %>%
      setNames(as.list(.), .)

    SusieOut <- purrr::map(selected_phenotypes, ~finemapPhenotype(., 
                                                                           selected_qtl_group, 
                                                                           GenotypeMatrix, 
                                                                           TrainCovariateSet, 
                                                                           CisDistance,
                                                                           variant_list = VariantList, 
                                                                           ))
    SusieRes <- purrr::map(SusieOut, extractResults) %>%
        purrr::transpose() %>% 
        CleanSusieData(RegionDf)
    HoldOutPredictions <- GetPredictions(SusieRes,
                                            GenotypeMatrix,
                                            TestGeneVector,
                                            TrainCovarModel,
                                            TestCovariateSet) 
    HoldOutPredictions
}





ResidualizeMolecularData <- function(GeneVector,
                                     Covariates,
                                     Coefs) {
SortedGeneVector <- GeneVector
rownames(SortedGeneVector) <- SortedGeneVector$sample_id
print(head(SortedGeneVector))
SortedGeneVector <- SortedGeneVector[rownames(Covariates),] 
TestResids <- data.frame(Observed = SortedGeneVector$phenotype_value - data.matrix(Covariates) %*% Coefs) %>% 
                tibble::rownames_to_column('sample_id') 
GeneVectorResidualized <- GeneVector %>% 
                    left_join(TestResids,by = 'sample_id')
GeneVectorResidualized
}




GetPredictions <- function(SusieRes,
                        GenotypeMatrix,
                        GeneVector,
                        TrainCoefs,
                        Covariates) {

Samples <- GeneVector %>% pull(sample_id)
VarNames <- data.frame(variant=rownames(GenotypeMatrix)) %>% 
            mutate(variant = str_replace(variant,'chrchr','chr')) %>%
            pull(variant)
rownames(GenotypeMatrix) <- VarNames
CleanedSusie <- SusieRes %>% 
        mutate(variant = str_replace(variant,'chrchr','chr')) 

GenotypeDataFinemappedVariants <- GenotypeMatrix[CleanedSusie$variant,Samples]
PredictedValues <- t(GenotypeDataFinemappedVariants) %*% as.vector(SusieRes %>% pull(posterior_mean)) %>% 
    data.frame() %>% 
    tibble::rownames_to_column('sample_id') %>% 
    dplyr::rename('Predicted' = 2)
ObservedValues <- GeneVector %>% ResidualizeMolecularData(Covariates,TrainCoefs)
MergedData <- PredictedValues %>% left_join(ObservedValues,by ='sample_id') 
MergedData
}


EstimateBetaHat <- function(MolecularVector,CovarMatrix) {
BetaHat <- solve(crossprod(CovarMatrix), crossprod(CovarMatrix,MolecularVector ))
BetaHat
}

GetFoldPCData <- function(FoldPCData,
                       Fold) {
FoldLabel <- paste0('Fold',Fold)
FoldData <- FoldPCData[[FoldLabel]]
FoldData
}

GetFoldTrainPCs <- function(FoldPCData,
                            Fold) {
FoldPCData <- GetFoldPCData(FoldPCData,Fold)
TrainPCs <- FoldPCData[['TrainPCs']]
TrainPCs
}

GetFoldTestPCs <- function(FoldPCData,
                            Fold) {
FoldPCData <- GetFoldPCData(FoldPCData,Fold)
TestPCs <- FoldPCData[['TestPCs']]
TestPCs
}

MergeMolecularGeneticPCs <- function(ExpressionPCs,GeneticPCs) {
GeneticPCsFiltered <- GeneticPCs %>%
        data.frame() %>% 
        tibble::rownames_to_column('ID') %>% 
        select(ID,contains('GENETIC'))
MergedData <- ExpressionPCs %>% 
                left_join(GeneticPCsFiltered,by = 'ID') %>% 
                tibble::column_to_rownames('ID')  %>% 
                data.matrix()
MergedData
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

CleanSusieData <- function(res,region_df) {
      
  #Extract information about all variants
  variant_df <- purrr::map_df(res$variant_df, identity, .id = "phenotype_id")
  if(nrow(variant_df) > 0){
    variant_df <- variant_df %>%
      dplyr::left_join(region_df, by = "phenotype_id") %>%
      tidyr::separate(variant_id, c("chr", "pos", "ref", "alt"),sep = "_", remove = FALSE) %>%
      dplyr::mutate(chr = stringr::str_remove_all(chr, "chr")) %>%
      dplyr::mutate(cs_index = cs_id) %>%
      dplyr::mutate(cs_id = paste(phenotype_id, cs_index, sep = "_"))
    pip_df = extractPipsFromVariantDf(variant_df)
    variant_df = dplyr::left_join(variant_df, pip_df, by = c("phenotype_id", "variant_id", "cs_id", "cs_index"))
  }
  
  #Extract lbf variable df and format correctly for export
  lbf_df_res <- purrr::map_df(res$lbf_df, identity, .id = "phenotype_id")
  if(nrow(lbf_df_res) > 0){
    lbf_df <- lbf_df_res %>%
      dplyr::left_join(region_df, by = "phenotype_id") %>%
      tidyr::separate(variant_id, c("chromosome", "position", "ref", "alt"),sep = "_", remove = FALSE) %>%
      dplyr::mutate(chromosome = stringr::str_remove_all(chromosome, "chr")) %>%
      dplyr::select(-ref, -alt) %>%
      dplyr::rename(molecular_trait_id = phenotype_id, variant = variant_id) %>%
      dplyr::select(molecular_trait_id, region, variant, chromosome, position, lbf_variable1:lbf_variable10)
    lbf_df <- lbf_df %>% dplyr::mutate(position = as.integer(position))
  } else {
    lbf_df = empty_lbf_df
  }
    cs_df <- purrr::map_df(res$cs_df, identity, .id = "phenotype_id")
    if(nrow(cs_df) > 0){
      cs_df = dplyr::left_join(cs_df, region_df, by = "phenotype_id") %>%
        dplyr::mutate(cs_index = cs_id) %>%
        dplyr::mutate(cs_id = paste(phenotype_id, cs_index, sep = "_")) %>%
        dplyr::transmute(molecular_trait_id = phenotype_id, cs_id, cs_index, region, cs_log10bf, cs_avg_r2, cs_min_r2, cs_size, low_purity)
      
      #Extract information about variants that belong to a credible set
      in_cs_variant_df <- dplyr::filter(variant_df, !is.na(cs_index) & !low_purity) %>%
        dplyr::transmute(molecular_trait_id = phenotype_id, variant = variant_id, chromosome = chr, position = pos, 
                         ref, alt, cs_id, cs_index, region, pip, z, cs_min_r2, cs_avg_r2, cs_size, posterior_mean, posterior_sd, cs_log10bf)
    } else{
      #Initialize empty tibbles with correct column names
      in_cs_variant_df = empty_in_cs_variant_df
      cs_df = empty_cs_df
    }


     

    #Extract information about all variants
    if(nrow(variant_df) > 0){
      variant_df_transmute <- dplyr::transmute(variant_df, molecular_trait_id = phenotype_id, variant = variant_id, 
              chromosome = chr, position = pos, ref, alt, cs_id, cs_index, low_purity, region, pip, z, posterior_mean, posterior_sd, X_column_scale_factors)  
      variant_df <- dplyr::bind_cols(variant_df_transmute, dplyr::select(variant_df,alpha1:mu2_10))
    } else{
      variant_df = empty_variant_df
    variant_df
    }
}
MergeCovars <- function(GeneticPCs,ExpressionPCs) {
Merged <- GeneticPCs %>% 
                tibble::rownames_to_column('ID') %>% 
                left_join(ExpressionPCs,by = 'ID') %>% 
                tibble::column_to_rownames('ID') %>% 
                data.matrix()
Merged
    
}



