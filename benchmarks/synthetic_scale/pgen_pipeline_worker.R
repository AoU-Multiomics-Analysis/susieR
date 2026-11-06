# Benchmark the current entry point with matched regional TSV and PGEN inputs.
# Keep all production calculations. Add only CLI arguments and timing hooks.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 6L, args[[3L]] %in% c("tabix_pipeline", "pgen_pipeline"))
input_dir <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = TRUE)
mode <- args[[3L]]
repo_dir <- normalizePath(args[[4L]], mustWork = TRUE)
expected_dimensions <- as.integer(args[5:6])
stopifnot(length(expected_dimensions) == 2L, !anyNA(expected_dimensions),
          all(expected_dimensions > 0L))
.libPaths(c(Sys.getenv("BENCH_LIBRARY"), .libPaths()))
suppressPackageStartupMessages(library(readr))

helper_dir <- file.path(output_dir, "helpers")
dir.create(helper_dir, showWarnings = FALSE)
helper_paths <- list.files(file.path(repo_dir, "R", "utils"),
                          pattern = "[.]R$", full.names = TRUE)
stopifnot(length(helper_paths) > 0L,
          all(file.copy(helper_paths, helper_dir, overwrite = TRUE)))
Sys.setenv(SUSIER_FUNCTIONS_PATH = helper_dir)

count_fixed <- function(text, pattern) {
  sum(gregexpr(pattern, text, fixed = TRUE)[[1L]] > 0L)
}
if (mode == "pgen_pipeline") {
  source(file.path(repo_dir, "benchmarks", "synthetic_scale", "pgen_region_reader.R"),
         local = .GlobalEnv)
  readSyntheticPgenDosage <- function(chr, start, end, dosage_file) {
    readPgenRegionHardcalls(
      pgen_path = dosage_file,
      pvar_path = file.path(dirname(dosage_file), "synthetic.pvar"),
      psam_path = file.path(dirname(dosage_file), "synthetic.psam"),
      # The synthetic cis window crosses the start of the chromosome.
      # PVAR positions are positive, so this selects the same Tabix records.
      chrom = chr, start = max(1, start), end = end, as_double = TRUE
    )
  }
  helper_path <- file.path(helper_dir, "SusieFunctions.R")
  helper <- paste(readLines(helper_path, warn = FALSE), collapse = "\n")
  original_reader <- "eQTLUtils::extractGenotypeMatrixFromDosage("
  stopifnot(count_fixed(helper, original_reader) == 2L)
  helper <- gsub(original_reader, "readSyntheticPgenDosage(", helper, fixed = TRUE)
  write_lines(helper, helper_path)
}

timing_rows <- list()
record_stage <- function(name, seconds) {
  timing_rows[[length(timing_rows) + 1L]] <<- tibble::tibble(
    stage = name, elapsed_seconds = seconds
  )
  write_tsv(dplyr::bind_rows(timing_rows), file.path(output_dir, "stages.tsv"))
}
genotype_path <- file.path(input_dir, if (mode == "pgen_pipeline") {
  "synthetic.pgen"
} else {
  "dosages.tsv.bgz"
})
stopifnot(file.exists(genotype_path))
pipeline_cli <- c(
  "--genotype_matrix", genotype_path,
  "--sample_meta", file.path(input_dir, "samples.tsv"),
  "--phenotype_list", file.path(input_dir, "permutations.tsv"),
  "--expression_matrix", file.path(input_dir, "phenotypes.bed"),
  "--covariates", file.path(input_dir, "covariates.tsv"),
  "--out_prefix", file.path(output_dir, "synthetic"),
  "--cisdistance", "1000000", "--MAF", "0"
)

script <- paste(readLines(file.path(repo_dir, "R", "scripts", "susie.R"),
                          warn = FALSE), collapse = "\n")
cli_line <- "opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))"
stopifnot(count_fixed(script, cli_line) == 1L)
script <- sub(cli_line, paste0(
  "opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list), ",
  "args = pipeline_cli)"
), script, fixed = TRUE)
reader_trace <- if (mode == "pgen_pipeline") {
  "trace('readSyntheticPgenDosage', where = .GlobalEnv,"
} else {
  "trace('extractGenotypeMatrixFromDosage', where = asNamespace('eQTLUtils'),"
}
hooks <- paste(c(
  "trace('LoadData', tracer = quote(load_started <- proc.time()[['elapsed']]),",
  "      exit = quote(record_stage('load_data', proc.time()[['elapsed']] - load_started)), print = FALSE)",
  reader_trace,
  "      tracer = quote(genotype_started <- proc.time()[['elapsed']]),",
  "      exit = quote(record_stage('genotype_read', proc.time()[['elapsed']] - genotype_started)), print = FALSE)",
  "trace('fitSusieModel', tracer = quote({",
  "        stopifnot(identical(dim(X), expected_dimensions))",
  "        saveRDS(.Random.seed, file.path(output_dir, 'rng_before.rds'))",
  "        fit_started <- proc.time()[['elapsed']]",
  "      }), exit = quote({",
  "        record_stage('fit', proc.time()[['elapsed']] - fit_started)",
  "        saveRDS(.Random.seed, file.path(output_dir, 'rng_after.rds'))",
  "      }), print = FALSE)"
), collapse = "\n")
source_line <- "source(file.path(FunctionPath, \"MemoryHelpers.R\"))"
stopifnot(count_fixed(script, source_line) == 1L)
script <- sub(source_line, paste(source_line, hooks, sep = "\n"), script, fixed = TRUE)
snapshot <- file.path(output_dir, "entrypoint.R")
write_lines(script, snapshot)
saveRDS(list(
  R = R.version.string,
  BLAS = extSoftVersion()[["BLAS"]],
  susieR = as.character(packageVersion("susieR")),
  thread_variables = Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")),
  mode = mode, genotype_file = basename(genotype_path)
), file.path(output_dir, "environment.rds"))
set.seed(20261106L)
elapsed <- system.time(source(snapshot, local = .GlobalEnv))[["elapsed"]]
record_stage("entrypoint_total", elapsed)
message(mode, " completed: ", elapsed, " s")
