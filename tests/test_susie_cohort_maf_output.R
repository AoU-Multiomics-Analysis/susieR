# Catch MAF attached by a repeated label, lost from cached output, or added to
# the full, LBF, or credible-set tables. The native fit supplies physical rows.
run_susie_cohort_maf_output_tests <- function(repo_dir = ".") {
  functions_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH",
                              unset = file.path(repo_dir, "R", "utils"))
  helpers <- new.env(parent = globalenv())
  for (file in c("InitFunctions.R", "SusieFunctions.R", "SusieOutputFunctions.R")) {
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
  n <- 120L
  p <- 20L
  variant_ids <- paste0("chr1_", 1000L + seq_len(p), "_A_G")
  variant_ids[[2L]] <- variant_ids[[1L]]
  X <- matrix(rnorm(n * p), nrow = n,
              dimnames = list(NULL, variant_ids))
  X[, 2L] <- X[, 1L] + rnorm(n, sd = 0.1)
  X[, 4L] <- X[, 3L] + rnorm(n, sd = 0.1)
  y <- cbind(0.65 * X[, 1L] - 0.55 * X[, 3L] + rnorm(n))
  fitted <- suppressMessages(susieR::susie(
    X, y, L = 10L, estimate_prior_variance = FALSE,
    scaled_prior_variance = 0.1, compute_univariate_zscore = TRUE
  ))
  fitted$variant_id <- variant_ids
  stopifnot("The fixture must have two native sets in physical predictor order" =
              identical(fitted$sets$cs, list(L1 = c(1L, 2L), L2 = c(3L, 4L))))

  # The extraction predates the metadata. Cached output must read the source fit.
  old_results <- list(gene_maf = fitted, gene_old = fitted)
  extractions <- purrr::map(old_results, helpers$extractResults)
  results <- old_results
  results$gene_maf$cohort_maf <- c(0.125, 0.375, NA_real_, 0.25, (5:20) / 100)
  regions <- tibble::tibble(phenotype_id = rev(names(results)), region = "chr1_1_2000")
  fits_before <- serialize(results, NULL)
  seed_before <- .Random.seed
  reference <- helpers$formatSusieOutputs(old_results, regions)
  tables <- helpers$formatSusieOutputs(results, regions)
  cached <- helpers$formatSusieOutputs(results, regions, extracted_results = rev(extractions))

  failures <- character()
  check <- function(condition, label) {
    if (!isTRUE(condition)) failures <<- c(failures, label)
  }
  compact <- tables$in_cs_variant_df
  expected_maf <- c(0.125, 0.375, NA_real_, 0.25, rep(NA_real_, 4L))
  check(nrow(compact) == 8L &&
          identical(compact$molecular_trait_id, rep(names(results), each = 4L)) &&
          identical(compact$variant, rep(variant_ids[1:4], 2L)),
        "Compact output must keep each physical CS row when variant labels repeat")
  check(is.double(compact$cohort_maf) && identical(compact$cohort_maf, expected_maf),
        "Compact MAF must use physical predictor positions and numeric NA for old fits")
  check(is.double(reference$in_cs_variant_df$cohort_maf) &&
          identical(reference$in_cs_variant_df$cohort_maf, rep(NA_real_, 8L)),
        "Compact output from all old fits must have numeric unknown MAF")
  check(identical(tables, cached),
        "Cached extraction must use MAF from its source fit in fitted-outcome order")
  check(identical(dplyr::select(compact, -dplyr::any_of("cohort_maf")),
                  dplyr::select(reference$in_cs_variant_df, -dplyr::any_of("cohort_maf"))),
        "MAF must preserve all existing compact values and rows")
  for (name in c("variant_df", "lbf_df", "cs_df", "covariate_explained_variants")) {
    check(!"cohort_maf" %in% names(tables[[name]]) &&
            identical(tables[[name]], reference[[name]]),
          paste(name, "must keep its existing columns, order, and values"))
  }
  full_columns <- c(
    "molecular_trait_id", "variant", "chromosome", "position", "ref", "alt",
    "cs_id", "cs_index", "low_purity", "region", "pip", "z", "posterior_mean",
    "posterior_sd", "X_column_scale_factors", paste0("alpha", 1:10),
    paste0("mu_", 1:10), paste0("mu2_", 1:10), "phenotype_residual_SD",
    "feature_converged"
  )
  check(identical(names(tables$variant_df), full_columns),
        "Full output must keep the exact alpha, mu, and mu2 column range")
  check(all(purrr::map_lgl(extractions, function(extracted) {
    all(purrr::map_lgl(extracted, ~ !"cohort_maf" %in% names(.x)))
  })), "Public extraction tables must keep their existing schemas")

  output_dir <- tempfile("susier-cohort-maf-output-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
  files <- helpers$writeSusieOutputs(results, tables, file.path(output_dir, "mixed"))
  output_tables <- c(summary = "in_cs_variant_df", full = "variant_df", lbf = "lbf_df")
  for (path_key in names(output_tables)) {
    name <- output_tables[[path_key]]
    check(identical(as.data.frame(arrow::read_parquet(files[[path_key]])),
                    as.data.frame(tables[[name]])),
          paste("Serialized", path_key, "must preserve every column and value"))
  }
  check(identical(arrow::read_parquet(files[["summary"]])$cohort_maf, expected_maf),
        "Compact Parquet must preserve numeric MAF values and missing values")
  check(identical(readRDS(files[["rds"]]), results) &&
          identical(serialize(results, NULL), fits_before),
        "Export must keep fit metadata and scientific values unchanged")
  check(identical(.Random.seed, seed_before), "Export must keep the RNG state unchanged")

  no_cs <- fitted
  no_cs$sets$cs <- list()
  no_cs$sets$cs_index <- integer()
  no_cs$sets$purity <- fitted$sets$purity[0, , drop = FALSE]
  empty_results <- list(list(), list(gene_maf = no_cs))
  empty_regions <- list(regions[0, ], regions[regions$phenotype_id == "gene_maf", ])
  for (i in seq_along(empty_results)) {
    empty <- helpers$formatSusieOutputs(empty_results[[i]], empty_regions[[i]])
    check(nrow(empty$in_cs_variant_df) == 0L &&
            identical(empty$in_cs_variant_df$cohort_maf, double()),
          "Empty compact output must include numeric MAF")
    check(all(purrr::map_lgl(empty[c("variant_df", "lbf_df", "cs_df")],
                            ~ !"cohort_maf" %in% names(.x))),
          "Other empty output tables must keep their existing schemas")
    empty_files <- helpers$writeSusieOutputs(empty_results[[i]], empty,
                                             file.path(output_dir, paste0("empty", i)))
    check(identical(arrow::read_parquet(empty_files[["summary"]])$cohort_maf, double()),
          "Empty compact Parquet must keep the numeric MAF type")
  }
  check(!"cohort_maf" %in% names(helpers$InitEmptyInCSVariantDf()),
        "The shared compact constructor default must keep the checkpointed schema")
  if (length(failures)) stop(paste(failures, collapse = "\n"), call. = FALSE)
  message("SuSiE compact cohort MAF tests passed: physical rows, cached output, old fits, empty schemas, and Parquet")
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
  run_susie_cohort_maf_output_tests(repo_dir)
}
