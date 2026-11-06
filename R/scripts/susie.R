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

# Run SuSiE fine-mapping for one scattered molecular-trait batch and emit
# credible-set, log-Bayes-factor, and full variant-level parquet outputs.

####### GET PATH TO FUNCTIONS ########
# The Docker image installs helper code in /opt/r/lib; local runs can override
# that path with SUSIER_FUNCTIONS_PATH.
FunctionPath <- Sys.getenv("SUSIER_FUNCTIONS_PATH", unset = "/opt/r/lib")

####### IMPORT FUNCTINS AND PARSE OTHER COMMAND LINE ARGUMENTS############ 
# Source shared data loaders, empty-output schemas, fine-mapping helpers, option
# definitions, and debug utilities.
source(file.path(FunctionPath, "ImportFunctions.R"))
source(file.path(FunctionPath, "InitFunctions.R"))
source(file.path(FunctionPath, "SusiePerformance.R"))
source(file.path(FunctionPath, "SusieFunctions.R"))
source(file.path(FunctionPath, "SusieOutputFunctions.R"))
source(file.path(FunctionPath, "OptParser.R"))
source(file.path(FunctionPath, "MemoryHelpers.R"))

message('Functions Loaded')



option_list <- Optlist()
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

# Load all WDL-provided data objects into the script environment using the
# shared loader so helper functions and entrypoint scripts agree on input names.
data_list <- LoadData(opt)
list2env(data_list,envir = environment())

########### INITIALIZE EMPTY DATAFRAMES #########
message('Init empty dataframes')
empty_variant_df <- InitEmptyVariantDf()
empty_lbf_df <- InitEmptyLbfDf()
empty_in_cs_variant_df <- InitEmptyInCSVariantDf(include_covariate_diagnostics = TRUE)
empty_cs_df <- InitEmptyCS()


######### LOAD DATA #######
message('summarizing  phenotype data')
# SummarizedExperiment keeps molecular traits, sample metadata, and assay values
# aligned while subsetting to the qtl_group used for fine-mapping.
se = eQTLUtils::makeSummarizedExperimentFromCountMatrix(assay = expression_matrix %>% 
                                                            select(-1,-2,-3) , 
                                                         row_data = phenotype_meta, 
                                                         col_data = sample_metadata, 
                                                         quant_method = "gene_counts",
                                                         reformat = FALSE)
selected_qtl_group = eQTLUtils::subsetSEByColumnValue(se, "qtl_group",'ALL')
################# RUN FINE MAPPING ##########

# Each Cromwell scatter currently passes a single batch, but the chunk helper is
# retained so the same code can scale if batching is reintroduced.
selected_chunk_group = splitIntoChunks(1, 1, length(unique(phenotype_list$group_id)))
selected_group_ids = unique(phenotype_list$group_id)[selected_chunk_group]

selected_phenotypes = phenotype_list %>%
  dplyr::filter(group_id %in% selected_group_ids) %>%
  dplyr::pull(phenotype_id) %>%
  setNames(as.list(.), .) 
variant_position_report_file <- paste0(opt$out_prefix, ".variant_position_summary.tsv")
if (file.exists(variant_position_report_file)) {
  file.remove(variant_position_report_file)
}
reuse_genotype_matrix <- isTRUE(opt$reuse_genotype_matrix)
reusable_genotype <- NULL
if (reuse_genotype_matrix && length(selected_phenotypes) > 1) {
  message("Reusable genotype matrix requested")
  reusable_genotype <- prepareReusableGenotype(selected_phenotypes,
                                               selected_qtl_group,
                                               genotype_file,
                                               covariates_matrix,
                                               cis_distance,
                                               ancestry_data = AncestryDf,
                                               MAF = MAF_threshold,
                                               variant_list = variant_list,
                                               additional_genotypes = additional_genotypes)
} else if (reuse_genotype_matrix) {
  message("Reusable genotype matrix requested for a single phenotype; using per-phenotype path")
}
message('Fine-mapping begin')
# Run finemapPhenotype for each selected trait and keep the raw result object as
# an RDS for debugging or post hoc inspection.
results = purrr::map(selected_phenotypes, ~finemapPhenotype(., selected_qtl_group, 
                                                              genotype_file, 
                                                              covariates_matrix, 
                                                              cis_distance,
                                                              ancestry_data = AncestryDf, 
                                                              MAF = MAF_threshold,
                                                              variant_list =variant_list,
                                                              additional_genotypes = additional_genotypes,
                                                              reusable_genotype = reusable_genotype,
                                                              variant_report_file = variant_position_report_file
                                                              ))
tables <- formatSusieOutputs(results, region_df)
output_files <- writeSusieOutputs(results, tables, opt$out_prefix)
