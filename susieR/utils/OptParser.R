Optlist <- function() {
option_list <- list(
  # TODO look around if there is a package recognizing delimiter in dataset
  optparse::make_option(c("--MAF"), default=NULL,
    help="Minor allele frequency filter applied to genotype matrix", metavar="type",type="numeric"),
  optparse::make_option(c("--AncestryMetadata"), type="character", default=NULL,
    help="File that contains ancestry grouping assignments to calculate MAF per population", metavar="type"),
  optparse::make_option(c("--phenotype_meta"), type="character", default=NULL,
    help="Phenotype metadata file path (tab separated)", metavar="type"),
  optparse::make_option(c("--sample_meta"), type="character", default=NULL,
    help="Sample metadata file path (tab separated)", metavar="type"),
  optparse::make_option(c("--cv_meta"), type="character", default=NULL,
    help="Sample metadata file path (tab separated)", metavar="type"),
  optparse::make_option(c("--expression_matrix"), type="character", default=NULL,
    help="Expression matrix file path (genes in rows, samples in columns)", metavar="type"),
  optparse::make_option(c("--phenotype_list"), type="character", default=NULL,
    help="Path to phenotype list file", metavar="type"),
  optparse::make_option(c("--genotype_matrix"), type="character", default=NULL,
    help="Genotype dosage matrix extracted from VCF", metavar="type"),
  optparse::make_option(c("--covariates"), type="character", default=NULL,
    help="Path to covariates file in QTLtools format", metavar="type"),
  optparse::make_option(c("--out_prefix"), type="character", default="./finemapping_output",
    help="Prefix of output files", metavar="type"),
  optparse::make_option(c("--qtl_group"), type="character", default=NULL,
    help="Value of the current qtl_group", metavar="type"),
  optparse::make_option(c("--cisdistance"), type="integer", default=1000000,
    help="Cis distance (bp) from center of gene [default \"%default\"]", metavar="number"),
  optparse::make_option(c("--chunk"), type="character", default="1 1",
    help="Chunking (e.g. '5 10' = 5th of 10 chunks) [default \"%default\"]", metavar="type"),
  optparse::make_option(c("--eqtlutils"), type="character", default=NULL,
    help="Optional path to eQTLUtils package [default \"%default\"]", metavar="type"),
  optparse::make_option(c("--write_full_susie"), type="character", default="true",
    help="If 'true' full SuSiE output will not be written. Set 'false' to write all. [default \"%default\"]", metavar="type"),
  optparse::make_option(c("--train_test_split"), type="numeric",
    help="float specifying how much of the data to use to train the models and how much to validate", metavar="type",default = NULL),
  optparse::make_option(c("--n_folds"), type="numeric",
    help="specifies number of folds to use", metavar="type",default = NULL),
  optparse::make_option(c("--VariantList"), type="character",
    help="file that contains gnomad common variants, will use this to filter genotype data if present", metavar="type",default = NULL)
  
)
option_list
}


