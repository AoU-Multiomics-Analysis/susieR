suppressPackageStartupMessages(library("devtools"))
suppressPackageStartupMessages(library("susieR"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("Rsamtools"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("arrow"))
if (!requireNamespace("Rfast", quietly = TRUE)) {
  message('Installing Rfast')
  install.packages("Rfast", repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages(library("Rfast"))

####### GET PATH TO FUNCTIONS ########
#FunctionsPathOptList <- list(
    #optparse::make_option(
                          #c("--FunctionsPath"),
                          #type='character',
                          #default = "/opt/r/lib"
    #)
#)
#FunctionPathOpt <- optparse::parse_args(optparse::OptionParser(option_list=FunctionsPathOptList))
#FunctionPath <- FunctionPathOpt$FunctionsPath

####### IMPORT FUNCTINS AND PARSE OTHER COMMAND LINE ARGUMENTS########## 
source(paste0(FunctionPath,'/ImportFunctions.R'))
source(paste0(FunctionPath,'/InitFunctions.R'))
source(paste0(FunctionPath,'/SusieCVFunctions.R'))
source(paste0(FunctionPath,'/OptParser.R'))
message('Functions Loaded')

option_list <- Optlist()
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

data_list <- LoadData(opt)
list2env(data_list,envir = environment())


########### INITIALIZE EMPTY DATAFRAMES #########

empty_variant_df <- InitEmptyVariantDf()
empty_lbf_df <- InitEmptyLbfDf()
empty_in_cs_variant_df <- InitEmptyInCSVariantDf()
empty_cs_df <- InitEmptyCS()



###### PARSE COMMAND LINE ARGUMENTS ########## 
data_list <- LoadData(opt)
list2env(data_list,envir = environment())

######## CHUNK DATA ##########
# convert sample list into sample metadata 
# required by eQTLUtils and assign sampels to train and holdout  
message('Loading sample metadata')
selected_chunk_group = splitIntoChunks(1, 1, length(unique(phenotype_list$group_id)))
selected_group_ids = unique(phenotype_list$group_id)[selected_chunk_group]
selected_phenotypes = phenotype_list %>%
  dplyr::filter(group_id %in% selected_group_ids) %>%
  dplyr::pull(phenotype_id) %>%
  setNames(as.list(.), .) 

########## BEGIN ############

sample_metadata_folds <- LoadSampleMetadataCV(opt$sample_meta,AncestryDf,nFolds =n_folds ) 

se = eQTLUtils::makeSummarizedExperimentFromCountMatrix(assay = expression_matrix %>% 
                                                                select(-1,-2,-3) , 
                                                             row_data = phenotype_meta, 
                                                             col_data = sample_metadata, 
                                                             quant_method = "gene_counts",
                                                             reformat = FALSE)
   
message('Extracting gene meta')
gene_meta = dplyr::filter(SummarizedExperiment::rowData(se) %>% as.data.frame(), phenotype_id == output_prefix)
gene_vector = eQTLUtils::extractPhentypeFromSE(output_prefix, se, "counts") %>%
    dplyr::mutate(phenotype_value_std = qnorm((rank(phenotype_value, na.last = "keep") - 0.5) / sum(!is.na(phenotype_value))))

#message('Fine-mapping begin')
genotype_matrix_full = eQTLUtils::extractGenotypeMatrixFromDosage(
    chr = gene_meta$chromosome, 
    start = gene_meta$phenotype_pos - cis_distance, 
    end = gene_meta$phenotype_pos + cis_distance, 
    dosage_file = genotype_file) 
genotype_type_matrix_one_percent <- genotype_matrix_full %>% 
        filterMAF(AncestryDf,MAF_threshold = 0.01,variant_list = variant_list)

######### RUN CROSS VALIDATION ANALYSIS AND FINE MAPPING ###########
R2_data <- data.frame()
for (k in c(1:n_folds)) {
    message(paste0('Running on fold:',k))
    sample_metadata <-  sample_metadata_folds %>%
            mutate(qtl_group = case_when(fold == k ~ 'holdout',TRUE ~ 'train'))
    holdout_samples <- sample_metadata %>% filter(qtl_group == 'holdout')
    train_sample <- sample_metadata %>% filter(qtl_group == 'train')
    subset_expression_matrix <- expression_matrix %>% 
                                select(1,2,3,4,all_of(train_sample$sample_id))
    PCA_for_fold <- compute_pcs(subset_expression_matrix) 

    covariate_matrix_subset <- covariates_matrix[train_sample$sample_id,] %>% 
        select(contains('GENETIC')) %>% 
        MergeCovars(PCA_for_fold)
    
    se = eQTLUtils::makeSummarizedExperimentFromCountMatrix(assay = subset_expression_matrix %>% 
                                                                select(-1,-2,-3) , 
                                                             row_data = phenotype_meta, 
                                                             col_data = sample_metadata, 
                                                             quant_method = "gene_counts",
                                                             reformat = FALSE)
    selected_phenotype <- output_prefix
    selected_qtl_group <- eQTLUtils::subsetSEByColumnValue(se, "qtl_group",'train')
    selected_phenotypes = phenotype_list %>%
      dplyr::filter(group_id %in% selected_group_ids) %>%
      dplyr::pull(phenotype_id) %>%
      setNames(as.list(.), .)
    FullSusie <- purrr::map(selected_phenotypes, ~finemapPhenotype(., selected_qtl_group, 
                                                                  genotype_matrix_full, 
                                                                  covariates_matrix, 
                                                                  cis_distance,
                                                                  AncestryDf,
                                                                  MAF = 0
                                                                  ))
    OnePercentAFSusie <- purrr::map(selected_phenotypes, ~finemapPhenotype(., 
                                                                           selected_qtl_group, 
                                                                           genotype_type_matrix_one_percent, 
                                                                           covariates_matrix, 
                                                                           cis_distance,
                                                                           AncestryDf,
                                                                           MAF = 0.01
                                                                           ))
    OnePercentRes <- purrr::map(OnePercentAFSusie, extractResults) %>%
        purrr::transpose() %>% 
        CleanSusieData(region_df)

    FullSusieRes <- purrr::map(FullSusie, extractResults) %>%
        purrr::transpose() %>% 
        CleanSusieData(region_df)


    hold_samples <- sample_metadata %>% filter(qtl_group == 'holdout')
    OnePercentHoldoutData <- GetPredictions(OnePercentRes,genotype_matrix_full,gene_vector,covariates_matrix) %>% 
                filter(sample_id %in% hold_samples$sample_id)
    FullSusieHoldoutData <- GetPredictions(FullSusieRes,genotype_matrix_full,gene_vector,covariates_matrix) %>% 
                filter(sample_id %in% hold_samples$sample_id)
    
    OnePercentSummary <- broom::glance(lm(Observed ~ Predicted,data =OnePercentHoldoutData)) %>% mutate(AF_threshold = 0.01)
    FullDataSummary <- broom::glance(lm(Observed ~ Predicted,data =FullSusieHoldoutData))%>% mutate(AF_threshold = 0)
    Out <- bind_rows(FullDataSummary,OnePercentSummary) %>% mutate(Fold = k,gene = opt$output_prefix)
    R2_data <- bind_rows(Out,R2_data) 
}

R2_data %>% write_tsv(outfile)
