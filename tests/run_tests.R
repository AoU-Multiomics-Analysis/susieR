# Run from the repository root, or pass the repository path as the first argument.
benchmark_library <- Sys.getenv("BENCH_LIBRARY", unset = "")
if (nzchar(benchmark_library)) {
  if (!dir.exists(benchmark_library)) {
    stop("BENCH_LIBRARY does not exist: ", benchmark_library, call. = FALSE)
  }
  .libPaths(c(benchmark_library, .libPaths()))
}

args <- commandArgs(trailingOnly = TRUE)
repo_dir <- normalizePath(if (length(args) > 0L) args[[1L]] else getwd(),
                          mustWork = TRUE)

required_packages <- c(
  "data.table", "dplyr", "tidyr", "stringr", "purrr", "readr", "susieR",
  "matrixStats", "Matrix", "Rsamtools", "eQTLUtils", "SummarizedExperiment", "arrow"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing test dependencies: ", paste(missing_packages, collapse = ", "),
       ". Use the SuSiE image, or set BENCH_LIBRARY for a local library.",
       call. = FALSE)
}

test_files <- c(
  "test_covariate_reader.R", "test_susie_performance.R", "test_susie_fit_caches.R",
  "test_marginal_stat_reuse.R", "test_susie_native_fallback.R", "test_finemap_paths.R",
  "test_covariate_explained_outputs.R", "test_susie_output_pip.R",
  "test_susie_phenotype_residual_sd.R",
  "test_susie_feature_convergence.R",
  "test_susie_cv_output_pip.R", "test_memory_debug.R"
)
for (test_file in test_files) {
  test_path <- file.path(repo_dir, "tests", test_file)
  if (!file.exists(test_path)) {
    stop("Missing test file: ", test_path, call. = FALSE)
  }
  source(test_path)
}

run_covariate_reader_tests(repo_dir)
run_susie_performance_tests(repo_dir)
run_susie_native_fallback_tests(repo_dir)
fit_helpers <- new.env(parent = globalenv())
functions_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH",
                             unset = file.path(repo_dir, "R", "utils"))
sys.source(file.path(functions_path, "SusiePerformance.R"), envir = fit_helpers)
if (isTRUE(fit_helpers$.susieFitSourceSupported(asNamespace("susieR")))) {
  run_susie_fit_cache_tests(repo_dir)
  run_marginal_stat_reuse_tests(repo_dir)
} else {
  message("Installed SuSiE source uses native fallback; cache-specific tests skipped.")
}
run_finemap_path_tests(repo_dir)
run_covariate_explained_output_tests(repo_dir)
run_susie_output_pip_tests(repo_dir)
run_susie_phenotype_residual_sd_tests(repo_dir)
run_susie_feature_convergence_tests(repo_dir)
if (file.exists(file.path(functions_path, "SusieCVFunctions.R"))) {
  run_susie_cv_output_pip_tests(repo_dir)
}
run_memory_debug_tests(repo_dir)
message("All SuSiE regression tests passed.")
