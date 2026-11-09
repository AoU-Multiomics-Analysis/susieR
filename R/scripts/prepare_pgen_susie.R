#!/usr/bin/env Rscript
benchmark_library <- Sys.getenv("BENCH_LIBRARY", unset = "")
if (nzchar(benchmark_library)) .libPaths(c(benchmark_library, .libPaths()))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
functions_dir <- Sys.getenv("SUSIER_FUNCTIONS_PATH", unset = "/opt/r/lib")
for (name in c("PgenFunctions.R", "PgenBatchFunctions.R", "PgenQcFunctions.R", "PgenPreparationFunctions.R")) {
  source(file.path(functions_dir, name))
}
option <- function(name, type = "character", default = NULL) {
  make_option(paste0("--", name), type = type, default = default,
    dest = gsub("-", "_", name, fixed = TRUE))
}
options <- list(option("pgen"), option("pvar"), option("psam"), option("worklists-file"),
  option("cis-distance", "integer", 1000000L), option("mac", "double", 0),
  option("min-call-rate", "double", .95), option("plink2", default = "plink2"),
  option("plink-memory-mb", "integer", 1024L), option("output-dir"))
opt <- parse_args(OptionParser(option_list = options))
required <- c("pgen", "pvar", "psam", "worklists_file", "output_dir")
missing <- required[vapply(required, function(name) is.null(opt[[name]]), logical(1))]
if (length(missing)) stop("Missing required arguments: ", paste(missing, collapse = ", "), call. = FALSE)
preparePgenBatches(opt)
