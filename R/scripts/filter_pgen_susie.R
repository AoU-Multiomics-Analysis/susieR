#!/usr/bin/env Rscript
benchmark_library <- Sys.getenv("BENCH_LIBRARY", unset = "")
if (nzchar(benchmark_library)) .libPaths(c(benchmark_library, .libPaths()))
suppressPackageStartupMessages(library(optparse))
function_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH", unset = "/opt/r/lib")
source(file.path(function_path, "PgenFunctions.R"))
source(file.path(function_path, "PgenQcFunctions.R"))
options <- list(
  make_option("--pgen", type = "character"),
  make_option("--pvar", type = "character"),
  make_option("--psam", type = "character"),
  make_option("--mac", type = "double", default = 0),
  make_option("--min-call-rate", type = "double", default = .95),
  make_option("--output-dir", type = "character"),
  make_option("--plink2", type = "character", default = "plink2"),
  make_option("--plink-memory-mb", type = "double", default = 1024)
)
opt <- parse_args(OptionParser(option_list = options), convert_hyphens_to_underscores = TRUE)
required <- c("pgen", "pvar", "psam", "output_dir")
missing <- required[vapply(required, function(name) is.null(opt[[name]]), logical(1))]
if (length(missing)) stop("Missing required arguments: ", paste(missing, collapse = ", "), call. = FALSE)
filterPgenFile(opt)
