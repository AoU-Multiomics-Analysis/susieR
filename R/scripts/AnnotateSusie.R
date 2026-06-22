library(tidyverse)
library(data.table)
library(rtracklayer)
library(arrow)
library(plyranges)
library(bedr)

# Annotate merged SuSiE fine-mapping output with gene, regulatory, conservation,
# and GVS/VAT variant annotations.

########## LOADING FUNCTIONS ##################

# Read gnomAD gene constraint metrics and convert pLI into a coarse constrained
# vs. unconstrained label for downstream joins by Ensembl gene ID.
load_constraint_data <- function(gnomad_constraint_data) {
message('Loading gnomad data')
message(paste0('Using file: ',basename(gnomad_constraint_data)))

gnomad_dat <- fread(gnomad_constraint_data) %>% 
            filter(canonical == TRUE) %>% 
            select(gene_id,`lof.pLI`) %>% 
            mutate(Constrained =  case_when(lof.pLI > 0.9 ~ 'Constrained',TRUE ~ 'Unconstrained')) 
gnomad_dat
}


# load encode hg38 regulatory element annotations
# Expand ENCODE cCRE records with multiple labels into one logical column per
# annotation type, then return them as GRanges for overlap joins.
load_ENCODE_data <- function(ENCODE_cCRES_path) {
message('Loading ENCODE data')
message(paste0('Using file: ',basename(ENCODE_cCRES_path)))
ENCODE_cres <- fread(ENCODE_cCRES_path) %>% 
    dplyr::rename('seqnames' =1,'start' =2,'end' =3,'type' =6)  %>% 
    separate_rows(type,sep = ',') %>% 
    mutate(value = TRUE) %>% 
	pivot_wider(names_from = type,values_from = value,values_fill = list(value = FALSE)) %>%   
	makeGRangesFromDataFrame(keep.extra = TRUE)

ENCODE_cres

}

# Read GVS/VAT variant annotations and create the same chromosome-position-ref-alt
# ID used in fine-mapping output.
load_gvs_VAT_data <- function(PathVAT) {
message('Loading VAT data')
message(paste0('Using file: ',basename(PathVAT)))
VariantAnnotations <- fread(PathVAT) %>% 
       mutate(ID = paste(chrom,pos,ref,alt,sep = '_'))
VariantAnnotations
}

# Convert GVS allele frequencies to minor allele frequencies and derive the
# highest-MAF subpopulation plus matching AC/AN counts for each variant.
convert_gvs_af_to_maf <- function(VariantAnnotations) {
message('Converting GVS allele frequencies to MAF')

gvs_subpops <- c('afr','amr','eas','eur','mid','sas','oth')
gvs_af_cols <- c('gvs_all_af',paste0('gvs_',gvs_subpops,'_af'),'gvs_max_af')
gvs_count_cols <- c(
    'gvs_all_ac','gvs_all_an','gvs_max_ac','gvs_max_an',
    paste0('gvs_',gvs_subpops,'_ac'),
    paste0('gvs_',gvs_subpops,'_an')
)

VariantAnnotations <- VariantAnnotations %>%
    mutate(across(any_of(c(gvs_af_cols,gvs_count_cols)), as.numeric)) %>%
    mutate(across(any_of(gvs_af_cols), ~case_when(
        is.na(.x) ~ NA_real_,
        .x > 0.5 ~ 1 - .x,
        TRUE ~ .x
    )))

subpop_af_cols <- paste0('gvs_',gvs_subpops,'_af')
present_subpop_af_cols <- subpop_af_cols[subpop_af_cols %in% names(VariantAnnotations)]

if (length(present_subpop_af_cols) == 0) {
    return(VariantAnnotations)
}

subpop_af_matrix <- as.matrix(VariantAnnotations[present_subpop_af_cols])
all_missing_af <- rowSums(!is.na(subpop_af_matrix)) == 0
subpop_af_matrix[is.na(subpop_af_matrix)] <- -Inf
max_af_index <- max.col(subpop_af_matrix,ties.method = 'first')
row_index <- seq_len(nrow(VariantAnnotations))
max_subpop <- str_remove_all(present_subpop_af_cols[max_af_index],'^gvs_|_af$')

VariantAnnotations$gvs_max_subpop <- ifelse(all_missing_af,NA_character_,max_subpop)
VariantAnnotations$gvs_max_af <- ifelse(all_missing_af,NA_real_,subpop_af_matrix[cbind(row_index,max_af_index)])
VariantAnnotations$gvs_max_ac <- NA_real_
VariantAnnotations$gvs_max_an <- NA_real_

for (subpop in unique(max_subpop[!all_missing_af])) {
    subpop_rows <- !all_missing_af & max_subpop == subpop
    ac_col <- paste0('gvs_',subpop,'_ac')
    an_col <- paste0('gvs_',subpop,'_an')

    if (ac_col %in% names(VariantAnnotations)) {
        VariantAnnotations$gvs_max_ac[subpop_rows] <- VariantAnnotations[[ac_col]][subpop_rows]
    }
    if (an_col %in% names(VariantAnnotations)) {
        VariantAnnotations$gvs_max_an[subpop_rows] <- VariantAnnotations[[an_col]][subpop_rows]
    }
}

VariantAnnotations
}


# load plink allele frequency data
# This loader is retained for optional PLINK AF annotation, although the current
# workflow uses GVS/VAT frequencies instead.
load_afreq_data <- function(afreq_path){
message('Loading allele frequency data')
message(paste0('Using file: ',basename(afreq_path)))

dat <- fread(afreq_path) %>% 
        dplyr::rename('variant' = 'ID') %>% 
        mutate(variant = str_replace(variant,':','_')) %>% 
        mutate(chrom = case_when(str_detect(variant,'chrchr') ~ str_remove(variant,'chr'),TRUE ~ variant)) %>% 
        dplyr::select(variant,ALT_FREQS)
 
dat     
    
}

# function to load in FANTOM5 data, note this data 
# as of right now is downloaded from the UCSC genome browser
# FANTOM5 promoter/enhancer intervals are converted to GRanges so they can be
# joined in the same overlap pipeline as ENCODE cCRE annotations.
load_FANTOM5_data <- function(FANTOM5_path) {
message('Loading FANTOM5 data')
message(paste0('Using file: ',basename(FANTOM5_path)))

FANTOM5_df <- fread(FANTOM5_path) %>% 
    dplyr::rename('seqnames' ='V1','start' = 'V2','end' ='V3') %>% 
    mutate(FANTOM5 = TRUE) %>% 
    select(seqnames:end,FANTOM5) %>% 
    makeGRangesFromDataFrame(keep.extra = TRUE)
FANTOM5_df
}


############### QUERY FUNCTIONS #################

# helper function to query tabix data from vep. 
# This is not the vep generated VCF but rather 
# a cleaned table that has typically been summarized 
# to the worst consequence
# Build tabix query ranges from fine-mapped variants and widen the returned VEP
# consequence labels into one logical column per consequence term.
query_tabix <- function(fm_res,tabix_path) {
message('Querying tabix data')

ranges <- fm_res %>%
    mutate(chromosome = str_remove_all(chromosome,'chr')) %>% 
    transmute(variant = paste0('chr',chromosome,':',as.numeric(position)-1,'-',position)) %>% 
    pull(variant)
tabix_res <- tabix(ranges,tabix_path,check.chr = FALSE) %>% 
        distinct()  %>%
        separate_rows(V5,sep='&') %>%  
        mutate(value = TRUE) %>% 
        mutate(V5 = str_remove_all(V5, "\\s+")) %>% 
        pivot_wider(names_from = V5,values_from = value,values_fill = FALSE) %>% 
        dplyr::rename('chromosome' = 'V1','start' = 'V2','ref' = 'V3','alt' = 'V4')
message('Tabix query sucessful')
tabix_res
}  

# helper function to query big wig file and 
# get scores
# Pull conservation scores at variant positions from a BigWig file and append
# them to the fine-mapping table.
query_bigwig <- function(fm_data,bw_path) {
message(paste0('Annotating with ',basename(bw_path)))


positions <- fm_data %>%
        mutate(start = position,end = position) %>%  	
	makeGRangesFromDataFrame()
bw <- BigWigFile(bw_path)
message('importing big wig')
phyloP_scores <- import(bw, selection = BigWigSelection(positions))

message('Import sucessful')

message('Merging bigwig data and variants')
hits <- findOverlaps(positions, phyloP_scores)
scores <- rep(NA_real_, length(positions))
scores[queryHits(hits)] <- mcols(phyloP_scores)$score[subjectHits(hits)]
fm_data$phylop <- scores
fm_data
}


# querys vep table and removes synonymous variant 
# annotations if another column is  true 
# Join preprocessed VEP consequences back to the fine-mapping table and suppress
# synonymous-only labels when a more severe consequence is also present.
query_vep_table <- function(fm_data,VEP_table) {
message('Annotating with VEP')
VEP_annotation_data <- query_tabix(fm_data ,VEP_table)

message('Merging VEP data with variants')
VEP_annotated <- fm_data %>% 
    mutate(start = as.character(position)) %>%
    mutate(chromosome = str_remove_all(chromosome,'chr')) %>% 
    mutate(chromosome = paste0('chr',chromosome)) %>%  
    left_join(VEP_annotation_data,by = c('chromosome','start','ref','alt')) %>% 
    mutate(
     synonymous_variant = ifelse(
      # check if ANY logical column other than keep is TRUE
      if_any(where(is.logical) & !matches("^synonymous_variant$"), ~ .x),
      FALSE,
      synonymous_variant
        )
      )
message('VEP merge sucessful')
VEP_annotated
}

# helper function to annotate fine mapped data with granges 
# object
# Convert variants into single-base GRanges and left-join overlapping annotation
# intervals, filling missing logical annotation flags with FALSE.
query_grange_data <- function(fm_data,grange_annotations) {
message('annotating data with grange')

annotated_fm_data <- fm_data %>% 
    mutate(start = position,end = position) %>% 
    mutate(chromosome = str_remove_all(chromosome,'chr')) %>% 
    mutate(chromosome = paste0('chr',chromosome)) %>% 
    makeGRangesFromDataFrame(keep.extra = TRUE) %>% 
    join_overlap_left(grange_annotations) %>% 
    data.frame() %>% 
    mutate(across(where(is.logical), ~replace_na(., FALSE))) %>% 
    dplyr::rename('chromosome' = 'seqnames')
annotated_fm_data

}



########### COMMAND LINE ARGUMENTS ########
message('Begin')
# Paths here are supplied by the WDL task. Optional sources that are not part of
# the current production annotation path remain commented for easy restoration.
option_list <- list(
  optparse::make_option(c("--SusieTSV"), type="character", default=NULL, metavar = "type"),
  optparse::make_option(c("--GencodeGTF"), type="character", default=NULL, metavar = "type"),
  optparse::make_option(c("--OutputPrefix"), type="character", default=NULL, metavar = "type"),
  #optparse::make_option(c("--PlinkAfreq"), type="character", default=NULL, metavar = "type"),
  optparse::make_option(c("--ENCODEcCRES"), type="character", default=NULL, metavar = "type"),
  #optparse::make_option(c("--VEPAnnotationsTable"), type="character", default=NULL, metavar = "type"),
  optparse::make_option(c("--gnomadConstraint"), type="character", default=NULL, metavar = "type"),
  optparse::make_option(c("--phyloPBigWig"), type="character", default=NULL, metavar = "type"),
  optparse::make_option(c("--FANTOM5"), type="character", default=NULL, metavar = "type"),
  optparse::make_option(c("--VAT"), type="character", default=NULL, metavar = "type")


)

message('Parsing command line arguments')
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

OutputPrefix <- opt$OutputPrefix


annotated_tsv <- paste0(opt$OutputPrefix,'_SusieMerged.annotated.tsv') 
message(paste0('Writing to:',annotated_tsv))

# this file needs to be a parsed version of the a VEP consequences VCF
#PathVEP <- opt$VEPAnnotationsTable
PathENCODE <- opt$ENCODEcCRES
PathPhyloP <- opt$phyloPBigWig
PathGnomad <- opt$gnomadConstraint
PathFANTOM5 <- opt$FANTOM5
PathSusie <- opt$SusieTSV
PathVAT <- opt$VAT
########### LOAD DATA ############

# Load all external annotation resources before touching the SuSiE data so
# missing reference files fail early.
FANTOM5_granges <- load_FANTOM5_data(PathFANTOM5)
ENCODE_data <- load_ENCODE_data(PathENCODE)
gnomad_data <- load_constraint_data(PathGnomad) 
VATData <- load_gvs_VAT_data(PathVAT)
#allele_frequencies <- load_afreq_data(opt$PlinkAfreq)

#susie_res <- load_finemapping_data(opt$SusieParquet)


message('Loading GTF')
# The Gencode GTF provides gene IDs, names, types, and strand-aware TSS
# positions used for distance-to-TSS annotations.
gene_data <- rtracklayer::readGFF(opt$GencodeGTF) %>% filter(type == 'gene')

message('Extracting TSS  locations')
tss_data <- gene_data %>% mutate(tss = case_when(strand == '+' ~ start,TRUE ~ end)) %>% 
            dplyr::select(seqid,tss,gene_id,gene_type,gene_name) %>% 
            mutate(start = tss,end = tss ) %>% 
            makeGRangesFromDataFrame(keep.extra = TRUE)



# annotate fine mapping data with internal data 
# generated from QTL mapping
message('Annotating fine-mapping data')
# Normalize variant IDs, attach gene/TSS metadata, and bin PIPs before adding
# external annotations.
annotated_fm_res <-  fread(PathSusie) %>%
  mutate(group = OutputPrefix) %>%
  mutate(variant = str_remove_all(variant,'chr')) %>% 
  mutate(variant = paste0('chr',variant)) %>% 
  #left_join(allele_frequencies,by = 'variant' ) %>% 
  #mutate(
	    #MAF = case_when(ALT_FREQS > .5 ~ 1 -ALT_FREQS,TRUE ~ ALT_FREQS),
  		#posterior_mean = case_when(ALT_FREQS > .5 ~ -posterior_mean,TRUE ~ posterior_mean),
 		#ref = case_when(ALT_FREQS > .5 ~ alt ,TRUE ~ ref),
		#alt = case_when(ALT_FREQS > .5 ~ ref ,TRUE ~ alt)
		#) %>% 
  #mutate(
  #      AF_bin = case_when(
  #        MAF  < 0.01 ~ "rare (0.1–1%)",
  #        MAF >= 0.01 & MAF < 0.05  ~ "low-freq (1–5%)",
  #        MAF >= 0.05 ~ "common (≥5%)"
  #      )
  #    ) %>% 
    mutate(gene_id = str_remove(molecular_trait_id,'.*_'))  %>% 
    left_join(tss_data %>% data.frame() %>% select(-seqnames,-start,-end,-width,-strand) ,by = 'gene_id')  %>% 
    mutate(distTSS = as.numeric(position) - as.numeric(tss),
           PIP_bin = cut(pip,breaks = 5)
    )  %>% 
    #dplyr::select(-ALT_FREQS ) %>%
    mutate(PIP_decile = cut(pip, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)) %>%  
    mutate(gene_id = str_remove(molecular_trait_id,'\\..*')) 

message('Annotating fine-mapping data with external data sources')
#drop_columns  <- c('alleles','locus','chrom','pos','ref.y','alt.y','AF','AC','AN','ALL_p_value_hwe','ALL_p_value_excess_het')
drop_columns  <- c('alleles','locus','chrom','pos','AF','AC','AN','ALL_p_value_hwe','ALL_p_value_excess_het','strand')

# Layer interval annotations, conservation scores, gene constraints, and VAT/GVS
# data into one final variant table.
full_annotated_data <- annotated_fm_res %>%
            query_grange_data(ENCODE_data) %>% 
            query_grange_data(FANTOM5_granges) %>% 
            query_bigwig(PathPhyloP) %>% 
            left_join(gnomad_data,by = 'gene_id') %>% 
            left_join(VATData,by = c('variant' = 'ID')) %>% 
            convert_gvs_af_to_maf() %>%
            dplyr::select(-drop_columns) %>% 
            dplyr::rename('ref' = 'ref.x','alt' = 'alt.x','ENCODE_ID_1' = 'V4','ENCODE_ID_2' = 'V5') 

message('Cleaning annotated data')
# Collapse overlapping promoter/enhancer signals into a single broad annotation
# label used by downstream summaries.
cleaned_full_annotated_data <- full_annotated_data %>% 
    #mutate(
    # synonymous_variant = ifelse(
      # check if ANY logical column other than keep is TRUE
     # if_any(where(is.logical) & !matches("^synonymous_variant$"), ~ .x),
     # FALSE,
     # synonymous_variant
     #   )
    #) %>% 
  mutate(Annotation = case_when(dELS == TRUE ~ 'Enhancer', pELS == TRUE ~ 'Enhancer',PLS == TRUE ~ 'Promoter',FANTOM5 == TRUE ~ 'Promoter')) 



message('Writing to output') 
cleaned_full_annotated_data %>% write_tsv(annotated_tsv)
