#!/usr/bin/env Rscript
benchmark_library <- Sys.getenv("BENCH_LIBRARY", unset = "")
if (nzchar(benchmark_library)) .libPaths(c(benchmark_library, .libPaths()))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
function_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH", unset = "/opt/r/lib")
source(file.path(function_path, "PgenBatchFunctions.R"))
options <- lapply(c("expected-outcomes", "summary-files", "credible-set-files", "index-files",
                    "diagnostic-files", "output-dir"), function(name) {
  make_option(paste0("--", name), type = "character")
})
opt <- parse_args(OptionParser(option_list = options), convert_hyphens_to_underscores = TRUE)
required <- c("expected_outcomes", "summary_files", "credible_set_files", "index_files", "diagnostic_files", "output_dir")
missing <- required[vapply(required, function(name) is.null(opt[[name]]), logical(1))]
if (length(missing)) stop("Missing required arguments: ", paste(missing, collapse = ", "), call. = FALSE)
read_file_list <- function(path) {
  paths <- readLines(pgenLocalFile(path, "Localized file list"), warn = FALSE)
  if (!length(paths) || anyNA(paths) || any(!nzchar(paths)) || any(grepl("[\r\n\t]", paths))) {
    stop("File lists need one nonempty localized path per line", call. = FALSE)
  }
  for (input in paths) pgenLocalFile(input, "Listed batch input")
  paths
}
message("Aggregating PGEN SuSiE batch summaries")
aggregatePgenBatches(
  expected_outcomes = opt$expected_outcomes,
  summary_files = read_file_list(opt$summary_files),
  credible_set_files = read_file_list(opt$credible_set_files),
  index_files = read_file_list(opt$index_files),
  diagnostic_files = read_file_list(opt$diagnostic_files),
  output_dir = opt$output_dir
)
message("PGEN batch aggregation complete")
