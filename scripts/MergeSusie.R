library(tidyverse)
library(data.table)
library(arrow)
library(optparse)

########### COMMAND LINE ARGUMENTS ########
option_list <- list(
  #TODO look around if there is a package recognizing delimiter in dataset
  optparse::make_option(c("--FilePaths"), type="character", default=NULL,
                        help="Phenotype metadata file path of genes used in expression-matrix. Tab separated", metavar = "type"),
  optparse::make_option(c("--OutputPrefix"), type="character", default=NULL,
                        help="Sample metadata file path of genes used in expression-matrix. Tab separated", metavar = "type")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

merged_parquet <- paste0(opt$OutputPrefix,'_SusieMerged.parquet') 
merged_tsv <- paste0(opt$OutputPrefix,'_SusieMerged.tsv.gz') 
aggregate_mode <- opt$SusieType

######## PARSE DATA #########
filepath_df <- fread(opt$FilePaths,header = FALSE) %>% dplyr::rename('path' = 1) %>% pull(path)
number_files <- filepath_df %>% length() 

message(paste0('Number of files found: ',number_files))
read_one_parquet <- function(path) {
  message("Reading: ", path)
  
  out <- tryCatch(
    arrow::read_parquet(path),
    error = function(e) {
      message("Skipping file due to read error: ", path)
      return(NULL)
    }
  )
  
  if (is.null(out)) {
    return(NULL)
  }
  
  out <- as.data.frame(out)
  
  if (nrow(out) == 0) {
    message("File is empty: ", path)
    return(NULL)
  }
  
  out
}

all_dat <- lapply(filepath_df, read_one_parquet)
all_dat <- Filter(Negate(is.null), all_dat)

if (length(all_dat) == 0) {
  stop("No non-empty parquet files were found.")
}

merged_df <- dplyr::bind_rows(all_dat)

data.table::fwrite(merged_df, merged_tsv)
arrow::write_parquet(merged_df, merged_parquet)

#lbf_df <- dplyr::tibble(
  #molecular_trait_id = character(),
  #region = character(),
  #variant = character(),
  #chromosome = character(),
  #position = integer(),
  #lbf_variable1 = numeric(),
  #lbf_variable2 = numeric(),
  #lbf_variable3 = numeric(),
  #lbf_variable4 = numeric(),
  #lbf_variable5 = numeric(),
  #lbf_variable6 = numeric(),
  #lbf_variable7 = numeric(),
  #lbf_variable8 = numeric(),
  #lbf_variable9 = numeric(),
  #lbf_variable10 = numeric()
#)

#variant_df <- dplyr::tibble(
  #molecular_trait_id = character(),
  #variant = character(),
  #chromosome = character(),
  #position = integer(),
  #ref = character(),
  #alt = character(),
  #cs_id = character(),
  #cs_index = character(),
  #region = character(),
  #pip = numeric(),
  #z = numeric(),
  #cs_min_r2 = numeric(),
  #cs_avg_r2 = numeric(),
  #cs_size = integer(),
  #posterior_mean = numeric(),
  #posterior_sd = numeric(),
  #cs_log10bf = numeric()
#)


#counter <- 0

#if (aggregate_mode == 'pip') { 
#for (x in filepath_df){
#current_dat <- arrow::read_parquet(x)
#variant_df <- bind_rows(variant_df,current_dat)
#counter <- counter + 1 
#print(counter)
    #} 
#}

#if (aggregate_mode == 'lbf') { 
#for (x in filepath_df){
#current_dat <- arrow::read_parquet(x)
#variant_df <- bind_rows(lbf_df,current_dat)
#counter <- counter + 1 
#print(counter)
    #} 
#}



#variant_df %>% fwrite(merged_tsv)
#variant_df %>% arrow::write_parquet(merged_parquet)
