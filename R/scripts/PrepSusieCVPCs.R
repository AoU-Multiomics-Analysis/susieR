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

# Prepare fold-specific expression PCs for SuSiE cross-validation. Inputs are
# injected by the WDL command block before this script is run.

FunctionPath <- Sys.getenv("SUSIER_FUNCTIONS_PATH", unset = "/opt/r/lib")
# Load the CV PC helpers and shared import utilities from the Docker image or a
# local override path.
source(file.path(FunctionPath, "SusieCVPrepCovars.R"))
source(file.path(FunctionPath, "ImportFunctions.R"))


######### LOAD DATA ################
# Build the sample/fold metadata used to estimate PCs within each training fold.
AncestryDf <- LoadAncestryData(AncestryPredictions)
SampleMetaDataFolds <- LoadSampleMetadataCV(SampleMetaDataPath,AncestryDf,nFolds = 10)
PositionTSS <- extract_TSS_pos(GTFPath)
SampleList <- fread(basename(SampleListPath),header = FALSE) %>% dplyr::rename('ID' = 1) %>% pull(ID)

# Restrict the expression matrix to localized samples and orient it as
# sample-by-feature before edgeR normalization.
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
# edgeR normalization is run before fold-level rank normalization in the helper
# functions.
DataCPM <- edgeR::cpm(DataEdgeR, log=FALSE) %>% data.frame() 

message('Normalizing CPMs')
# Return to BED-like feature ordering for the downstream PC helper.
NormalizedCPMs <- DataCPM %>% 
                t() %>% 
                data.frame() %>% 
                arrange(seqnames,start) %>%  
                dplyr::rename('#chr' = 'seqnames')

######## COMPUTE PCS ACROSS FOLDS ########
# Compute train/test PCs separately by fold to avoid leakage from held-out
# samples into the learned PC basis.
FoldPCData <- ComputePCsCV(SampleMetaDataFolds,ExpressionDf)
