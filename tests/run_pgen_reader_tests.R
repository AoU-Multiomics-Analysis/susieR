# Run from the repository root, or pass its path as the first argument.
benchmark_library <- Sys.getenv("BENCH_LIBRARY", unset = "")
if (nzchar(benchmark_library)) {
  if (!dir.exists(benchmark_library)) {
    stop("BENCH_LIBRARY does not exist: ", benchmark_library, call. = FALSE)
  }
  .libPaths(c(benchmark_library, .libPaths()))
}
args <- commandArgs(trailingOnly = TRUE)
repo_dir <- normalizePath(if (length(args)) args[[1L]] else getwd(),
                          mustWork = TRUE)
required_packages <- c("pgenlibr", "readr", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing test dependencies: ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}
source(file.path(repo_dir, "tests", "test_pgen_reader.R"))
run_pgen_reader_tests(repo_dir)
source(file.path(repo_dir, "tests", "test_pgen_native_imputation.R"))
run_pgen_native_imputation_tests(repo_dir)
message("All PGEN reader tests passed.")
