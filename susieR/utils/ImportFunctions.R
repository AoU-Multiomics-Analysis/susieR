ImportVariantList <- function(VariantListPath) {
    VariantList <- readr::read_tsv(VariantListPath) %>% 
                        dplyr::rename('variants'  = 1) %>% 
                        pull(variants)
    VariantList
}

importQtlmapCovariates <- function(covariates_path,load_colnames = FALSE){
  pc_matrix = read.table(covariates_path, check.names = F, header = T, stringsAsFactors = F)
  col_names <-  pc_matrix %>% pull(ID)
  pc_transpose = t(pc_matrix[,-1])
  colnames(pc_transpose) = pc_matrix$SampleID
  pc_df = dplyr::mutate(as.data.frame(pc_transpose), genotype_id = rownames(pc_transpose)) %>%
    dplyr::as_tibble() %>% 
    dplyr::select(genotype_id, dplyr::everything())
  
  #Make PCA matrix
  pc_matrix = as.matrix(dplyr::select(pc_df,-genotype_id))
  rownames(pc_matrix) = pc_df$genotype_id
  if (load_colnames == TRUE) {
      colnames(pc_matrix) <- col_names
  }
  return(pc_matrix)
}
importQtlmapPermutedPvalues <- function(perm_path){
  tbl = read.table(perm_path, check.names = F, header = T, stringsAsFactors = F) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(p_fdr = p.adjust(pval_beta, method = "fdr")) %>%
    dplyr::mutate(group_id = phenotype_id)
  return(tbl)
}


# Groups samples by ancestry and creates folds for cross validation
CreateCVFoldsMetadata <- function(SampleMetadataPath,AncestryDf,nFolds = 5) {
require(caret)
SampleMetadata <- readr::read_tsv(SampleMetadataPath) %>% 
        dplyr::rename('sample_id' =1 ) %>% 
        mutate(genotype_id = sample_id,qtl_group = 'ALL') %>% 
        mutate(sample_id = as.numeric(sample_id),genotype_id = as.character(genotype_id))  %>%
        left_join(AncestryDf,by = c('sample_id' = 'research_id')) %>% 
        group_by(ancestry_pred_other) %>%
        group_modify(~ {
            .x$fold <- createFolds(seq_len(nrow(.x)), k = nFolds, list = FALSE)
            .x
        }) %>%
        ungroup()
SampleMetadata
}

LoadAncestryData <- function(AncestryPath) {
message('Loading ancestry data')
AncestryDat <- readr::read_tsv(AncestryPath) %>% 
    select(research_id,ancestry_pred_other) %>%
    mutate(research_id = as.numeric(research_id))
AncestryDat
}


LoadCVData <- function(CVobj) {
#require(qs2)
message('Loading CV PCs and metadata ')
CVDat <- readRDS(CVobj)
}

LoadData <- function(opt_list) {
    require(data.table)
    message('Loading molecular data')
    
    if (is.null(opt_list$expression_matrix)) {
        stop('molecular data file is missing')
    }
    if (is.null(opt_list$phenotype_list)) {
        stop('tensorQTL file is missing')
    }
    if (is.null(opt_list$covariates)) {
        stop('Covariates file is missing')
    }
    if (is.null(opt_list$genotype_matrix)) {
        stop('Genotype file is missing')
    }
    
 
    expression_matrix = fread(opt_list$expression_matrix) %>% dplyr::rename('phenotype_id' = 'gene_id')

    message('Loading covariates')
    
    if(!is.null(opt_list$cv_meta)) {
        covariates_matrix = importQtlmapCovariates(opt_list$covariates,load_colnames = TRUE)
    } else {
        covariates_matrix = importQtlmapCovariates(opt_list$covariates,load_colnames = FALSE)
    }
    exclude_cov = apply(covariates_matrix, 2, sd) != 0
    covariates_matrix = covariates_matrix[,exclude_cov]

    # convert bed file into phenotype metadata required by eQTLUtils
    phenotype_meta<- expression_matrix %>% 
        select(1,2,3,4) %>% 
        dplyr::rename('chromosome' = 1,'phenotype_pos' = 2) %>% 
        mutate(strand  = 1) %>% 
        mutate(gene_id = phenotype_id,group_id = phenotype_id) %>% 
        select(phenotype_id,group_id,gene_id,chromosome,phenotype_pos,strand)
    # convert sample list into sample metadata 
    # required by eQTLUtils 
    message('Loading sample metadata')
    if(is.null(opt_list$cv_meta) && is.null(opt_list$sample_meta)) {
          stop('Sample metadata missing')
        }  
    if (!is.null(opt_list$sample_meta)) {
            sample_metadata <-  readr::read_tsv(opt_list$sample_meta) %>% 
                dplyr::rename('sample_id' =1 ) %>% 
                mutate(genotype_id = sample_id,qtl_group = 'ALL') %>% 
                mutate(sample_id = as.character(sample_id),genotype_id = as.character(genotype_id))
            }else {
            sample_metadata <- NULL
        }
    if (!is.null(opt_list$cv_meta)) {
            cv_meta <- LoadCVData(opt_list$cv_meta)  
                        
            }else {
            cv_meta <- NULL
        }
        # import permutation p values from tensorQTL. Note that i had 
    # to make slight changes to this function for it to work 
    message('Loading QTL stats')
    
    
    phenotype_table = importQtlmapPermutedPvalues(opt_list$phenotype_list)

    filtered_list = dplyr::filter(phenotype_table, p_fdr < 0.05) 
    phenotype_list = dplyr::semi_join(data.frame(group_id=phenotype_table$phenotype_id,phenotype_id=phenotype_table$phenotype_id) , filtered_list, by = "group_id")
    message("Number of phenotypes included for analysis: ", nrow(phenotype_list))
    #Keep only those phenotypes that are present in the expression matrix
    phenotype_list = dplyr::filter(phenotype_list, phenotype_id %in% expression_matrix$phenotype_id)

    if (!is.null(opt_list$AncestryMetadata)) {
        AncestryDf <- LoadAncestryData(opt_list$AncestryMetadata)
    }else {
        AncestryDf <- NULL
    }
    #if (!is.null(opt_list$n_folds)) {
        #n_folds <- as.numeric(opt_list$n_folds) 
        #}else {
        #n_folds <- NULL
    #}
   
    genotype_file <- opt_list$genotype_matrix
    cis_distance <- as.numeric(opt_list$cisdistance)
    output_prefix <- opt_list$out_prefix
    n_folds <- opt_list$n_folds
    variant_list <- opt_list$VariantList
    MAF_threshold <- opt_list$MAF
    region_df <- phenotype_meta %>% 
                    filter(phenotype_id %in% phenotype_list$phenotype_id) %>% 
                    transmute(phenotype_id,region = paste0(chromosome,':',
                                                           phenotype_pos - cis_distance,'-',
                                                           phenotype_pos + cis_distance))
    OutList <- list(
        covariates_matrix = covariates_matrix,
        phenotype_table = phenotype_table,
        phenotype_meta = phenotype_meta,
        filtered_list = filtered_list,
        phenotype_list = phenotype_list,
        expression_matrix = expression_matrix,
        variant_list = variant_list,
        MAF_threshold = MAF_threshold,
        AncestryDf = AncestryDf,
        output_prefix = output_prefix,
        cis_distance = cis_distance,
        genotype_file = genotype_file,
        sample_metadata = sample_metadata,
        region_df = region_df,
        cv_meta = cv_meta
        )
    OutList
}

