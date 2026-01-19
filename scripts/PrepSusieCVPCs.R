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


source('/opt/r/lib/SusieCVPrepCovars.R')
source('/opt/r/lib/ImportFunctions.R')


######### LOAD DATA ################
AncestryDf <- LoadAncestryData(AncestryPredictions)
SampleMetaDataFolds <- LoadSampleMetadataCV(SampleMetaDataPath,AncestryDf,nFolds = 10)
ExpressionDf <- fread(ExpressionBed)


######## COMPUTE PCS ACROSS FOLDS ########
FoldPCData <- ComputePCsCV(SampleMetaDataFolds,ExpressionDf)
