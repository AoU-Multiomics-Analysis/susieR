#!/usr/bin/env Rscript
benchmark_library <- Sys.getenv("BENCH_LIBRARY", unset = "")
if (nzchar(benchmark_library)) .libPaths(c(benchmark_library, .libPaths()))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
function_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH", unset = "/opt/r/lib")
source(file.path(function_path, "ImportFunctions.R"))
source(file.path(function_path, "PgenBatchFunctions.R"))
source(file.path(function_path, "PgenQcFunctions.R"))
options <- list(
  make_option("--phenotype-bed", type = "character"),
  make_option("--qtl-results", type = "character"),
  make_option("--chromosomes-file", type = "character"),
  make_option("--pgen-manifest", type = "character"),
  make_option("--cis-distance", type = "double"),
  make_option("--output-dir", type = "character"),
  make_option("--batch-size", type = "double", default = 50),
  make_option("--mac", type = "double", default = 0),
  make_option("--min-call-rate", type = "double", default = 0.95),
  make_option("--gene-map", type = "character", default = NULL),
  make_option("--select-top-phenotype-per-cluster", action = "store_true", default = FALSE),
  make_option("--top-phenotype-pvalue-column", type = "character", default = "qval")
)
opt <- parse_args(OptionParser(option_list = options), convert_hyphens_to_underscores = TRUE)
required <- c("phenotype_bed", "qtl_results", "cis_distance", "output_dir")
missing <- required[vapply(required, function(name) is.null(opt[[name]]), logical(1))]
if (length(missing)) stop("Missing required arguments: ", paste(missing, collapse = ", "), call. = FALSE)
validatePgenQcThresholds(opt$mac, opt$min_call_rate)
if (is.null(opt$pgen_manifest) == is.null(opt$chromosomes_file)) {
  stop("Supply exactly one of --pgen-manifest or --chromosomes-file", call. = FALSE)
}
manifest <- if (!is.null(opt$pgen_manifest)) readPgenManifest(opt$pgen_manifest) else NULL
chromosomes <- if (is.null(manifest)) {
  readLines(pgenLocalFile(opt$chromosomes_file, "Chromosome list"), warn = FALSE)
} else manifest$chromosome
message("Planning PGEN SuSiE batches")
plan <- planPgenBatches(
  phenotype_bed = opt$phenotype_bed, qtl_results = opt$qtl_results,
  chromosomes = chromosomes, cis_distance = opt$cis_distance, output_dir = opt$output_dir,
  batch_size = opt$batch_size, gene_map = opt$gene_map,
  select_top_phenotype_per_cluster = opt$select_top_phenotype_per_cluster,
  top_phenotype_pvalue_column = opt$top_phenotype_pvalue_column
)
if (!is.null(manifest)) {
  writePgenManifest(manifest, file.path(opt$output_dir, "chromosome_pgens.tsv"))
  writePgenQcPlan(manifest, plan$chromosomes, opt$output_dir)
}
message("PGEN batch planning complete")
