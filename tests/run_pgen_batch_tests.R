# Run from the repository root, or pass its path as the first argument.
benchmark_library <- Sys.getenv("BENCH_LIBRARY", unset = "")
if (nzchar(benchmark_library)) .libPaths(c(benchmark_library, .libPaths()))
args <- commandArgs(trailingOnly = TRUE)
repo_dir <- normalizePath(if (length(args)) args[[1L]] else getwd(), mustWork = TRUE)
source(file.path(repo_dir, "tests", "test_pgen_batches.R"))
run_pgen_batch_tests(repo_dir)
message("All PGEN batch tests passed.")
