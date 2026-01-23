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
PositionTSS <- extract_TSS_pos(GTFPath)
SampleList <- fread(basename(SampleListPath),header = FALSE) %>% dplyr::rename('ID' = 1) %>% pull(ID)


CountsData <- fread(ExpressionBed)
CountDataTransposed <- CountData %>%
    dplyr::select(-Description) %>% 
    column_to_rownames('Name') %>%
    dplyr::select(any_of(SampleList)) %>% 
    t() %>% 
    data.frame()

DataEdgeR <- edgeR::DGEList(CountDataFiltered)
DataEdgeR <- edgeR::calcNormFactors(DataEdgeR)

message('Computing CPMs')
DataCPM <- edgeR::cpm(DataEdgeR, log=FALSE) %>% data.frame() 

message('Normalizing CPMs')
NormalizedCPMs <- DataCPM %>% 
                t() %>% 
                data.frame() %>% 
                arrange(seqnames,start) %>%  
                dplyr::rename('#chr' = 'seqnames')

######## COMPUTE PCS ACROSS FOLDS ########
FoldPCData <- ComputePCsCV(SampleMetaDataFolds,ExpressionDf)
