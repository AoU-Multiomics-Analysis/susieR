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


######## PARSE DATA #########
filepath_df <- fread(opt$)
