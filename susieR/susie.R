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
FunctionsPathOptList <- list(
    optparse::make_option(
                          c("--FunctionsPath"),
                          type='character',
                          default = "/opt/r/lib"
    )
)
FunctionPathOpt <- optparse::parse_args(optparse::OptionParser(option_list=FunctionsPathOptList))
FunctionPath <- FunctionPathOpt$FunctionsPath

####### IMPORT FUNCTINS AND PARSE OTHER COMMAND LINE ARGUMENTS############ 
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


######### LOAD DATA #######
se = eQTLUtils::makeSummarizedExperimentFromCountMatrix(assay = expression_matrix %>% 
                                                            select(-1,-2,-3) , 
                                                         row_data = phenotype_meta, 
                                                         col_data = sample_metadata, 
                                                         quant_method = "gene_counts",
                                                         reformat = FALSE)
selected_qtl_group = eQTLUtils::subsetSEByColumnValue(se, "qtl_group",'ALL')
################# RUN FINE MAPPING ##########

selected_chunk_group = splitIntoChunks(1, 1, length(unique(phenotype_list$group_id)))
selected_group_ids = unique(phenotype_list$group_id)[selected_chunk_group]

selected_phenotypes = phenotype_list %>%
  dplyr::filter(group_id %in% selected_group_ids) %>%
  dplyr::pull(phenotype_id) %>%
  setNames(as.list(.), .) 
message('Fine-mapping begin')
results = purrr::map(selected_phenotypes, ~finemapPhenotype(., selected_qtl_group, 
                                                              genotype_file, 
                                                              covariates_matrix, 
                                                              cis_distance,
                                                              MAF = MAF_threshold,
                                                              variant_list = VariantList
                                                              ))

#Only proceed if the there are more than 0 phenotypes
message("Number of overall unique group_ids: ", length(unique(phenotype_list$group_id)))
message("Number of groups in the batch: ", length(selected_group_ids))
message("Number of phenotypes in the batch: ", length(selected_phenotypes))

#Extract credible sets from finemapping results
message(" # Extract credible sets from finemapping results")
res = purrr::map(results, extractResults) %>%
    purrr::transpose()
  
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
  
#Extract information about credible sets
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
}

if (nrow(variant_df) == 0 && nrow(cs_df) == 0 && nrow(in_cs_variant_df) == 0) {
  arrow::write_parquet(in_cs_variant_df, paste0(opt$out_prefix, ".parquet"))
  arrow::write_parquet(empty_lbf_df, paste0(opt$out_prefix, ".lbf_variable.parquet"))
  arrow::write_parquet(empty_variant_df, paste0(opt$out_prefix, ".full_susie.parquet"))
  message("There are no credible sets. Write empty matrices and stop execution.")
  quit(save = "no", status = 0)
} 

in_cs_variant_df <- in_cs_variant_df %>% dplyr::mutate(position = as.integer(position))

# find how many unique phenotypes there are per gene
in_cs_variant_gene_df <- in_cs_variant_df %>% 
  dplyr::left_join(phenotype_meta %>% dplyr::select(phenotype_id, gene_id, group_id), by = c("molecular_trait_id" = "phenotype_id")) %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::mutate(uniq_phenotypes_count = length(base::unique(molecular_trait_id))) %>% 
  dplyr::ungroup()

## if it is gene expression write full sumstats
#if (all(in_cs_variant_gene_df$molecular_trait_id == in_cs_variant_gene_df$gene_id) | opt$write_full_susie) {
arrow::write_parquet(in_cs_variant_df, paste0(opt$out_prefix, ".parquet"))
arrow::write_parquet(lbf_df, paste0(opt$out_prefix, ".lbf_variable.parquet"))
arrow::write_parquet(variant_df, paste0(opt$out_prefix, ".full_susie.parquet"))











#} else { 
  ## generate connected components per gene
#message("Building connected components!")
#susie_cc <- make_connected_components_from_cs(susie_all_df = in_cs_variant_gene_df, cs_size_threshold = 200)
#needed_phenotype_ids <- susie_cc$molecular_trait_id %>% base::unique()
#  arrow::write_parquet(in_cs_variant_df, paste0(opt$out_prefix, ".parquet"))
#  arrow::write_parquet(lbf_df, paste0(opt$out_prefix, ".lbf_variable.parquet"))
#  arrow::write_parquet(variant_df, paste0(opt$out_prefix, ".full_susie.parquet"))




#in_cs_variant_df_filt <- in_cs_variant_df %>% dplyr::filter(molecular_trait_id %in% needed_phenotype_ids)
#cs_df_filt <- cs_df %>% dplyr::filter(molecular_trait_id %in% needed_phenotype_ids)
#variant_df_filt <- variant_df %>% dplyr::filter(molecular_trait_id %in% needed_phenotype_ids)
#lbf_df_filt <- lbf_df %>% dplyr::filter(molecular_trait_id %in% needed_phenotype_ids)
#arrow::write_parquet(in_cs_variant_df_filt, paste0(opt$out_prefix, ".parquet"))
#arrow::write_parquet(lbf_df_filt, paste0(opt$out_prefix, ".lbf_variable.parquet"))
#arrow::write_parquet(variant_df_filt, paste0(opt$out_prefix, ".full_susie.parquet"))


