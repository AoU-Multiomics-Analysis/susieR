# Catch a scale taken from fitted noise variance or lost during output caching.
run_susie_phenotype_residual_sd_tests <- function(repo_dir = ".") {
  functions_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH",
                              unset = file.path(repo_dir, "R", "utils"))
  helpers <- new.env(parent = globalenv())
  for (file in c("SusiePerformance.R", "InitFunctions.R", "SusieFunctions.R",
                 "SusieOutputFunctions.R")) {
    sys.source(file.path(functions_path, file), envir = helpers)
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(20261109L)
  n <- 80L
  p <- 12L
  variant_ids <- paste0("chr1_", 1000L + seq_len(p), "_A_G")
  X <- matrix(rnorm(n * p), nrow = n,
              dimnames = list(NULL, variant_ids))
  # This final input vector has mean 7 and centered sum of squares 400.
  # Its SD differs from both the population SD and the fitted noise SD.
  y <- cbind(7 + rep(c(-3, -1, 1, 3), 20L))
  X[, 1L] <- drop(y) + rnorm(n, sd = 0.5)
  expected_sd <- sqrt(400 / 79)
  capture_fit <- function(fit) {
    set.seed(383L)
    invisible(capture.output(value <- suppressMessages(fit())))
    list(fit = value, rng = .Random.seed)
  }
  same_sd <- function(value) {
    is.double(value) && length(value) > 0L &&
      all(is.finite(value)) && all(abs(value - expected_sd) < 1e-14)
  }
  native_fit <- function(legacy = FALSE) {
    fitted <- susieR::susie(X, y, L = 10,
      estimate_residual_variance = TRUE, estimate_prior_variance = TRUE,
      scaled_prior_variance = 0.1, verbose = TRUE,
      compute_univariate_zscore = legacy, min_abs_corr = 0.5)
    if (!legacy) fitted$z <- helpers$computeMarginalZScores(X, y)
    fitted$covariate_explained_variants <- helpers$emptyCovariateExplainedVariants()
    fitted
  }
  fits <- list()
  for (legacy in c(FALSE, TRUE)) {
    actual <- capture_fit(function() helpers$fitSusieModel(X, y, legacy = legacy))
    stopifnot("The common fit boundary must save the sample SD of its final y" =
                same_sd(actual$fit$phenotype_residual_SD))
    reference <- capture_fit(function() native_fit(legacy))
    without_scale <- actual$fit
    without_scale$phenotype_residual_SD <- NULL
    stopifnot("Recording phenotype SD must preserve the complete native fit" =
                identical(without_scale, reference$fit),
              "Recording phenotype SD must preserve RNG state" =
                identical(actual$rng, reference$rng),
              "The fixture must distinguish phenotype SD from fitted noise SD" =
                abs(sqrt(actual$fit$sigma2) - expected_sd) > 0.5)
    fits[[if (legacy) "legacy" else "individual"]] <- actual$fit
  }

  output_dir <- tempfile("susier-phenotype-sd-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
  regions <- tibble::tibble(phenotype_id = "gene1", region = "chr1_1_2000")
  for (name in names(fits)) {
    fitted <- fits[[name]]
    fitted$variant_id <- variant_ids
    stopifnot("The fixture must yield a credible set" = length(fitted$sets$cs) > 0L)
    extracted <- helpers$extractResults(fitted)
    stopifnot("Extraction must repeat the recorded SD for every fitted variant" =
                nrow(extracted$variant_df) == p &&
                same_sd(extracted$variant_df$phenotype_residual_SD))
    results <- list(gene1 = fitted)
    tables <- helpers$formatSusieOutputs(results, regions)
    cached <- helpers$formatSusieOutputs(results, regions,
      extracted_results = list(gene1 = extracted))
    stopifnot("Cached output must retain the phenotype SD" = identical(tables, cached))
    files <- helpers$writeSusieOutputs(results, tables, file.path(output_dir, name))
    saved <- readRDS(files[["rds"]])
    stopifnot("The saved fit must retain the phenotype SD" =
                same_sd(saved$gene1$phenotype_residual_SD))
    for (table in c("summary", "full")) {
      reloaded <- arrow::read_parquet(files[[table]])
      expected_tail <- c("phenotype_residual_SD",
                         if (table == "summary") "cohort_maf", "feature_converged")
      stopifnot("Both Parquets must repeat a numeric phenotype SD" =
                  nrow(reloaded) > 0L && same_sd(reloaded$phenotype_residual_SD),
                "Recorded phenotype scale must precede the convergence annotation" =
                  identical(tail(names(reloaded), length(expected_tail)), expected_tail))
    }
  }

  old_fitted <- fits$individual
  old_fitted$variant_id <- variant_ids
  old_fitted$phenotype_residual_SD <- NULL
  old_fitted$sigma2 <- 123
  old_tables <- helpers$formatSusieOutputs(list(gene1 = old_fitted), regions)
  old_files <- helpers$writeSusieOutputs(list(gene1 = old_fitted), old_tables,
                                        file.path(output_dir, "old"))
  for (table in c("summary", "full")) {
    values <- arrow::read_parquet(old_files[[table]])$phenotype_residual_SD
    stopifnot("Old saved fits must use numeric NA without inferring from sigma2" =
                is.double(values) && length(values) > 0L && all(is.na(values)))
  }
  empty_tables <- helpers$formatSusieOutputs(list(), regions[0, ])
  empty_files <- helpers$writeSusieOutputs(list(), empty_tables,
                                          file.path(output_dir, "empty"))
  for (table in c("summary", "full")) {
    values <- arrow::read_parquet(empty_files[[table]])$phenotype_residual_SD
    stopifnot("Empty Parquets must retain a numeric phenotype SD column" =
                identical(values, numeric()))
  }
  message("Phenotype residual SD tests passed: fit paths, outputs, cached extraction, old fits, and empty schemas")
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  benchmark_library <- Sys.getenv("BENCH_LIBRARY", unset = "")
  if (nzchar(benchmark_library)) .libPaths(c(benchmark_library, .libPaths()))
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(purrr)
  })
  args <- commandArgs(trailingOnly = TRUE)
  repo_dir <- normalizePath(if (length(args)) args[[1L]] else getwd(), mustWork = TRUE)
  run_susie_phenotype_residual_sd_tests(repo_dir)
}
