#!/usr/bin/env Rscript
benchmark_library <- Sys.getenv("BENCH_LIBRARY", unset = "")
if (nzchar(benchmark_library)) .libPaths(c(benchmark_library, .libPaths()))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
function_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH", unset = "/opt/r/lib")
source(file.path(function_path, "PgenBatchFunctions.R"))
metadata_names <- c("gene_rds_uris_file", "gene_parquet_uris_file", "gene_tsv_uris_file",
                     "gene_log_uris_file", "gene_provenance_uris_file")
pgen_metadata_names <- c("pgen_chromosomes_file", "pgen_uris_file", "pvar_uris_file", "psam_uris_file")
options <- lapply(c("expected-outcomes", "summary-files", "credible-set-files", "index-files",
                    "diagnostic-files", "full-susie-files", "outcome-index", "output-dir",
                    gsub("_", "-", c(metadata_names, pgen_metadata_names))), function(name) {
  make_option(paste0("--", name), type = "character")
})
options <- c(options, list(make_option("--full-only", action = "store_true", default = FALSE,
  help = "Combine existing full-SuSiE files using --outcome-index, without compact aggregation")))
opt <- parse_args(OptionParser(option_list = options), convert_hyphens_to_underscores = TRUE)
read_file_list <- function(path) {
  paths <- readLines(pgenLocalFile(path, "Localized file list"), warn = FALSE)
  if (!length(paths) || anyNA(paths) || any(!nzchar(paths)) || any(grepl("[\r\n\t]", paths))) {
    stop("File lists need one nonempty localized path per line", call. = FALSE)
  }
  for (input in paths) pgenLocalFile(input, "Listed batch input")
  paths
}
if (isTRUE(opt$full_only)) {
  incompatible <- c("expected_outcomes", "summary_files", "credible_set_files", "index_files",
    "diagnostic_files", metadata_names, pgen_metadata_names)
  if (any(!vapply(opt[incompatible], is.null, logical(1)))) {
    stop("--full-only cannot combine compact aggregation inputs or output URI metadata", call. = FALSE)
  }
  required <- c("outcome_index", "full_susie_files", "output_dir")
  missing <- required[vapply(required, function(name) is.null(opt[[name]]), logical(1))]
  if (length(missing)) stop("Missing required full-only arguments: ", paste(missing, collapse = ", "), call. = FALSE)
  message("Aggregating existing full-SuSiE files")
  aggregatePgenFullSusie(opt$outcome_index, read_file_list(opt$full_susie_files),
    file.path(opt$output_dir, "full_susie.parquet"))
  message("Full-SuSiE aggregation complete")
  quit(status = 0L)
}
if (!is.null(opt$outcome_index)) stop("--outcome-index requires --full-only", call. = FALSE)
required <- c("expected_outcomes", "summary_files", "credible_set_files", "index_files", "diagnostic_files", "output_dir")
missing <- required[vapply(required, function(name) is.null(opt[[name]]), logical(1))]
if (length(missing)) stop("Missing required arguments: ", paste(missing, collapse = ", "), call. = FALSE)
has_metadata <- !vapply(opt[metadata_names], is.null, logical(1))
if (any(has_metadata) && !all(has_metadata)) {
  stop("Supply all five gene-output URI metadata lists together", call. = FALSE)
}
has_pgen_metadata <- !vapply(opt[pgen_metadata_names], is.null, logical(1))
if (any(has_pgen_metadata) && !all(has_pgen_metadata)) {
  stop("Supply all four PGEN output metadata lists together", call. = FALSE)
}
read_uri_metadata <- function(path) {
  # Open the localized list file only. Its contents describe completed outputs;
  # they must never be treated as files to open or localize in this task.
  values <- readLines(pgenLocalFile(path, "Output URI metadata list"), warn = FALSE)
  if (!length(values) || anyNA(values) || any(!nzchar(values)) || any(grepl("[\r\n\t]", values))) {
    stop("Output URI metadata lists need one nonempty location per line", call. = FALSE)
  }
  values
}
metadata <- if (all(has_metadata)) {
  setNames(lapply(opt[metadata_names], read_uri_metadata), sub("_file$", "", metadata_names))
} else NULL
pgen_metadata <- if (all(has_pgen_metadata)) {
  setNames(lapply(opt[pgen_metadata_names], read_uri_metadata),
           c("chromosomes", "pgen_uris", "pvar_uris", "psam_uris"))
} else NULL
if (!is.null(pgen_metadata)) {
  expected <- pgenReadTextTable(opt$expected_outcomes)
  pgenRequireColumns(expected, c("batch_id", "chromosome"), "Expected outcomes")
  expected_batches <- expected |>
    dplyr::distinct(batch_id, chromosome)
  if (any(vapply(pgen_metadata, length, integer(1)) != nrow(expected_batches))) {
    stop("Each PGEN output metadata list length must match the expected batch count", call. = FALSE)
  }
  if (!identical(pgenAutosome(pgen_metadata$chromosomes), pgenAutosome(expected_batches$chromosome))) {
    stop("PGEN output metadata chromosomes do not match the expected batch order", call. = FALSE)
  }
}
message("Aggregating PGEN SuSiE batch summaries")
files <- aggregatePgenBatches(
  expected_outcomes = opt$expected_outcomes,
  summary_files = read_file_list(opt$summary_files),
  credible_set_files = read_file_list(opt$credible_set_files),
  index_files = read_file_list(opt$index_files),
  diagnostic_files = read_file_list(opt$diagnostic_files),
  output_dir = opt$output_dir
)
if (!is.null(opt$full_susie_files)) {
  message("Aggregating localized full-SuSiE gene tables")
  files[["full_susie"]] <- aggregatePgenFullSusie(files[["index"]],
    read_file_list(opt$full_susie_files), file.path(opt$output_dir, "full_susie.parquet"))
}
if (!is.null(metadata)) {
  message("Writing chromosome gene-output manifests from completed output metadata")
  do.call(writeChromosomeGeneManifests, c(list(
    outcome_index_file = files[["index"]],
    output_dir = file.path(opt$output_dir, "gene_manifests")), metadata))
}
if (!is.null(pgen_metadata)) {
  message("Writing a reusable PGEN manifest from completed output metadata")
  do.call(writePgenOutputManifest, c(pgen_metadata,
    list(output_file = file.path(opt$output_dir, "pgen_manifest.tsv"))))
}
message("PGEN batch aggregation complete")
