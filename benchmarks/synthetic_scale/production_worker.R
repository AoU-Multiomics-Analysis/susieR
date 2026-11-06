# Run snapshots of the actual updated entry point and helpers, without replacing
# production calculations. Only CLI injection and timing hooks are added.
.libPaths(c(Sys.getenv("BENCH_LIBRARY"), .libPaths()))
suppressPackageStartupMessages(library(readr))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 6L, args[[3]] == "production")
input_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- normalizePath(args[[2]], mustWork = TRUE)
repo_dir <- normalizePath(args[[4]], mustWork = TRUE)
expected_dimensions <- as.integer(args[5:6])
helper_dir <- file.path(output_dir, "helpers")
dir.create(helper_dir, showWarnings = FALSE)
stopifnot(all(file.copy(list.files(file.path(repo_dir, "R", "utils"), pattern = "[.]R$", full.names = TRUE),
                         helper_dir, overwrite = TRUE)))
Sys.setenv(SUSIER_FUNCTIONS_PATH = helper_dir)
timing_rows <- list()
record_stage <- function(name, seconds) {
  timing_rows[[length(timing_rows) + 1L]] <<- data.frame(stage = name, elapsed_seconds = seconds)
  write_tsv(dplyr::bind_rows(timing_rows), file.path(output_dir, "stages.tsv"))
}
production_cli <- c("--genotype_matrix", file.path(input_dir, "dosages.tsv.bgz"),
                    "--sample_meta", file.path(input_dir, "samples.tsv"),
                    "--phenotype_list", file.path(input_dir, "permutations.tsv"),
                    "--expression_matrix", file.path(input_dir, "phenotypes.bed"),
                    "--covariates", file.path(input_dir, "covariates.tsv"),
                    "--out_prefix", file.path(output_dir, "synthetic"),
                    "--cisdistance", "1000000", "--MAF", "0")
script <- paste(readLines(file.path(repo_dir, "R", "scripts", "susie.R")), collapse = "\n")
old <- "opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))"
stopifnot(sum(gregexpr(old, script, fixed = TRUE)[[1L]] > 0L) == 1L)
script <- sub(old, "opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list), args = production_cli)",
              script, fixed = TRUE)
hook <- paste(c(
  "trace('LoadData', tracer = quote(load_started <- proc.time()[['elapsed']]),",
  "      exit = quote(record_stage('load_data', proc.time()[['elapsed']] - load_started)), print = FALSE)",
  "trace('fitSusieModel', tracer = quote({",
  "        stopifnot(identical(dim(X), expected_dimensions))",
  "        saveRDS(.Random.seed, file.path(output_dir, 'rng_before.rds'))",
  "        fit_started <- proc.time()[['elapsed']]",
  "      }), exit = quote({",
  "        record_stage('fit', proc.time()[['elapsed']] - fit_started)",
  "        saveRDS(.Random.seed, file.path(output_dir, 'rng_after.rds'))",
  "      }), print = FALSE)"), collapse = "\n")
source_line <- "source(file.path(FunctionPath, \"MemoryHelpers.R\"))"
stopifnot(sum(gregexpr(source_line, script, fixed = TRUE)[[1L]] > 0L) == 1L)
script <- sub(source_line, paste(source_line, hook, sep = "\n"), script, fixed = TRUE)
snapshot <- file.path(output_dir, "entrypoint.R")
write_lines(script, snapshot)
saveRDS(list(R = R.version.string, BLAS = extSoftVersion()[["BLAS"]],
             susieR = as.character(packageVersion("susieR")),
             thread_variables = Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS"))),
        file.path(output_dir, "environment.rds"))
set.seed(20261106L)
profile_path <- Sys.getenv("SUSIER_R_PROFILE", unset = "")
if (nzchar(profile_path)) {
  Rprof(profile_path, interval = 0.02, memory.profiling = TRUE, gc.profiling = TRUE)
}
elapsed <- system.time(source(snapshot, local = .GlobalEnv))[["elapsed"]]
if (nzchar(profile_path)) Rprof(NULL)
record_stage("entrypoint_total", elapsed)
message("Actual production entry point completed: ", elapsed, " s")
