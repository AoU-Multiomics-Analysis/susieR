ComputePCs <- function(expression_df){

subsetted_expression_dat <- expression_df %>% select(-c(1,2,3,4))
pca_standardized <- PCAtools::pca(subsetted_expression_dat)
n_pcs <- PCAtools::chooseGavishDonoho( subsetted_expression_dat ,  var.explained = pca_standardized$sdev^2, noise = 1)
message(paste0('Using' , n_pcs,' PCs'))
pca_standardized$PCA_cleaned <- pca_standardized$rotated %>% 
   data.frame() %>%
   select(1:n_pcs) %>% 
   rownames_to_column('ID') %>% 
   mutate(ID = str_remove(ID,'X'))
pca_standardized$n_pcs <- n_pcs
pca_standardized
}

ComputeTrainPCs <- function(BedData,SampleMetadataCV,Fold) {
message('Computing Training data PCs fold:',Fold)
TrainMetadata <- SampleMetadataCV %>% filter(fold != Fold)
TrainSamples <- TrainMetadata %>% mutate(sample_id=as.character(sample_id)) %>% pull(sample_id)
TrainBedData <- BedData %>% 
                    select(1,2,3,4,all_of(as.character(TrainSamples)))
TrainPCs <- ComputePCs(TrainBedData)
TrainPCs
}

ComputeTestPCs <- function(BedData,
                            TrainPCs,
                           SampleMetadataCV,
                           Fold) {
message('Computing Test data PCs fold:',Fold)
nPCsTrain <- TrainPCs$n_pcs
TestMetaData <- SampleMetadataCV %>% filter(fold == Fold)
TestSamples <- TestMetaData %>% mutate(sample_id=as.character(sample_id)) %>% pull(sample_id)
TestBedData <- BedData %>% 
                    select(all_of(as.character(TestSamples)))
TestPCs <- t(data.matrix(TestBedData)) %*% data.matrix(TrainPCs$loadings) %>% data.frame() %>%select(1:nPCsTrain)
TestPCs
}

ComputePCsCV <- function(SampleMetadata,ExpressionDf) {
nFolds <- max(SampleMetadata$fold)
FoldLabels <- paste0("Fold",seq_len(nFolds))
FoldData <- c()
for (k in 1:nFolds) {
message('Running fold:',k)
TrainPCs <- ComputeTrainPCs(ExpressionDf,
                                  SampleMetadata,
                                 k)
TestPCs <- ComputeTestPCs(ExpressionDf,
                                  TrainPCs,
                                  SampleMetadata,
                                  k)  
FoldPCs <- list(TrainPCs = TrainPCs,TestPCs =TestPCs,fold = k )  
FoldData <- c(FoldPCs,FoldData)
    }
names(FoldData) <- FoldLabels
FoldData
}
