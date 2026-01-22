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
#suppressPackageStartupMessages(library("qs2"))

FunctionPath <- '/opt/r/lib'
####### IMPORT FUNCTINS AND PARSE OTHER COMMAND LINE ARGUMENTS########## 
source(paste0(FunctionPath,'/ImportFunctions.R'))
source(paste0(FunctionPath,'/InitFunctions.R'))
source(paste0(FunctionPath,'/SusieCVFunctions.R'))
source(paste0(FunctionPath,'/SusieFunctions.R'))
source(paste0(FunctionPath,'/OptParser.R'))
message('Functions Loaded')

option_list <- Optlist()
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

data_list <- LoadData(opt)
list2env(data_list,envir = environment())

if (is.null(cv_meta)) {
    stop('Missing CV data missing')
}

# Load Sample metadata from CV PC object
SampleMetaData <- cv_meta[['Metadata']] %>% 
            #dplyr::rename('sample_id' =1 ) %>% 
            mutate(genotype_id = sample_id,qtl_group = 'ALL') %>% 
            mutate(sample_id = as.character(sample_id),genotype_id = as.character(genotype_id))
nFolds <- max(as.numeric(SampleMetaData$fold))

OutputFile <- paste0(output_prefix,'_SusiePredictions.tsv')

########### INITIALIZE EMPTY DATAFRAMES #########

empty_variant_df <- InitEmptyVariantDf()
empty_lbf_df <- InitEmptyLbfDf()
empty_in_cs_variant_df <- InitEmptyInCSVariantDf()
empty_cs_df <- InitEmptyCS()


######## CHUNK DATA ##########
# convert sample list into sample metadata 
# required by eQTLUtils and assign sampels to train and holdout  
selected_chunk_group = splitIntoChunks(1, 1, length(unique(phenotype_list$group_id)))
selected_group_ids = unique(phenotype_list$group_id)[selected_chunk_group]
selected_phenotypes = phenotype_list %>%
  dplyr::filter(group_id %in% selected_group_ids) %>%
  dplyr::pull(phenotype_id) %>%
  setNames(as.list(.), .) 

########## FORMAT DATA AND LAOD GENOTYPE MATRIX ############

#sample_metadata_folds <- LoadSampleMetadataCV(opt$sample_meta,AncestryDf,nFolds =n_folds ) 

se = eQTLUtils::makeSummarizedExperimentFromCountMatrix(assay = expression_matrix %>% 
                                                                select(-1,-2,-3) , 
                                                             row_data = phenotype_meta, 
                                                             col_data = SampleMetaData, 
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

######### RUN CROSS VALIDATION ANALYSIS AND FINE MAPPING ###########

Predictions <- data.frame()
message('Running CV on all variants')
for (k in c(1:nFolds)) {
    message(paste0('Running on fold:',k))
    FoldPredicitons <- RunFoldCV(SampleMetaData,
                                        covariates_matrix,
                                        cv_meta,
                                        expression_matrix,
                                        genotype_matrix_full,
                                        output_prefix
                                        phenotype_meta,
                                        cis_distance,
                                        region_df,
                                        k
                                        ) %>%
                                        mutate(AF_threshold = 0,Fold = k)
    Predictions <- bind_rows(FoldPredicitons,Predictions)    
}

message('Filtering Genotype matrix')
genotype_type_matrix_one_percent <- genotype_matrix_full %>% 
        filterMAF(AncestryDf,variant_list = variant_list)
rm(genotype_matrix_full)


message('Running CV on 1% variants')
for (k in c(1:nFolds)) {
    message(paste0('Running on fold:',k))
    FoldPredicitons <- RunFoldCV(SampleMetaData,
                                        covariates_matrix,
                                        cv_meta,
                                        expression_matrix,
                                        genotype_matrix_full,
                                        output_prefix
                                        phenotype_meta,
                                        variant_list,
                                        cis_distance,
                                        region_df,
                                        k
                                        ) %>%
                                        mutate(AF_threshold = 0.01,Fold = k)
    Predictions <- bind_rows(FoldPredicitons,Predictions)    
}


Predictions %>% mutate(Gene = output_prefix) %>% write_tsv(OutputFile)

#OnePercentSummary <- broom::glance(lm(Observed ~ Predicted,data =OnePercentHoldoutData)) %>% mutate(AF_threshold = 0.01)
#FullDataSummary <- broom::glance(lm(Observed ~ Predicted,data =FullSusieHoldoutData))%>% mutate(AF_threshold = 0)
#Out <- bind_rows(FullDataSummary,OnePercentSummary) %>% mutate(Fold = k,gene = output_prefix)
#R2_data <- bind_rows(Out,R2_data) 


