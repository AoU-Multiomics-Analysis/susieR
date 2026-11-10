# Catch convergence flags taken from another outcome, purity, or output presence.
run_susie_feature_convergence_tests <- function(repo_dir = ".") {
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

  # This native fit has multiple credible sets and unassigned variants. Copying
  # only its convergence metadata keeps the scientific values fixed.
  set.seed(20261109L)
  n <- 120L
  p <- 20L
  variant_ids <- paste0("chr1_", 1000L + seq_len(p), "_A_G")
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
  fitted$phenotype_residual_SD <- stats::sd(drop(y))
  stopifnot("The fixture must converge and produce multiple credible sets" =
              isTRUE(fitted$converged) && length(fitted$sets$cs) >= 2L,
            "The fixture must retain unassigned variants" =
              length(setdiff(seq_len(p), unlist(fitted$sets$cs))) > 0L)

  results <- list(gene_false_first = fitted, gene_true = fitted,
                  gene_false_last = fitted, gene_old = fitted)
  all_converged <- results
  results$gene_false_first$converged <- FALSE
  results$gene_false_last$converged <- FALSE
  results$gene_old$converged <- NULL
  expected <- tibble::tibble(
    molecular_trait_id = names(results),
    feature_converged = c(FALSE, TRUE, FALSE, NA)
  )
  regions <- tibble::tibble(phenotype_id = rev(names(results)), region = "chr1_1_2000")
  tables <- helpers$formatSusieOutputs(results, regions)
  reference <- helpers$formatSusieOutputs(all_converged, regions)
  extracted <- purrr::map(results, helpers$extractResults)
  stopifnot("Cached extractions must not need convergence metadata" =
              all(purrr::map_lgl(extracted, function(outcome) {
                !any(purrr::map_lgl(outcome, ~ "feature_converged" %in% names(.x)))
              })))
  cached <- helpers$formatSusieOutputs(results, regions, extracted_results = rev(extracted))

  failures <- character()
  check <- function(condition, label) {
    if (!isTRUE(condition)) failures <<- c(failures, label)
  }
  same_flags <- function(table) {
    index <- match(table$molecular_trait_id, expected$molecular_trait_id)
    is.logical(table$feature_converged) && !anyNA(index) &&
      identical(table$feature_converged, expected$feature_converged[index])
  }
  table_names <- c("cs_df", "in_cs_variant_df", "variant_df")
  for (name in table_names) {
    table <- tables[[name]]
    check(nrow(table) > 0L && same_flags(table),
          paste(name, "must repeat the native convergence value for its own feature"))
    check(identical(tail(names(table), 1L), "feature_converged"),
          paste(name, "must append feature_converged after all existing fields"))
    check(identical(dplyr::select(table, -dplyr::any_of("feature_converged")),
                    dplyr::select(reference[[name]], -dplyr::any_of("feature_converged"))),
          paste(name, "must retain every row and scientific value when convergence changes"))
  }
  check(nrow(tables$cs_df) == length(results) * length(fitted$sets$cs),
        "Nonconverged features must retain all their credible sets")
  check(nrow(tables$variant_df) == length(results) * p,
        "All features must retain every fitted variant")
  for (name in c("in_cs_variant_df", "variant_df")) {
    table <- tables[[name]]
    index <- match(table$variant, variant_ids)
    check(!anyNA(index) && identical(unname(table$pip), unname(fitted$pip[index])),
          paste(name, "must retain native variant PIPs"))
    check(identical(table$phenotype_residual_SD,
                    rep(fitted$phenotype_residual_SD, nrow(table))),
          paste(name, "must retain phenotype_residual_SD"))
  }
  check(identical(tables, cached),
        "Cached extractions without convergence fields must match fresh output")
  check(identical(tables$lbf_df, reference$lbf_df) &&
          identical(tables$covariate_explained_variants,
                    reference$covariate_explained_variants),
        "Convergence metadata must preserve the remaining output tables")

  output_dir <- tempfile("susier-feature-convergence-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
  files <- helpers$writeSusieOutputs(results, tables, file.path(output_dir, "mixed"))
  check(identical(readRDS(files[["rds"]]), results),
        "The saved native fits must preserve convergence metadata and credible sets")
  for (name in c("summary", "full")) {
    table <- arrow::read_parquet(files[[name]])
    check(nrow(table) > 0L && same_flags(table),
          paste("Serialized", name, "must preserve logical convergence values and missing values"))
    check(identical(tail(names(table), 1L), "feature_converged"),
          paste("Serialized", name, "must append feature_converged"))
  }

  empty_tables <- helpers$formatSusieOutputs(list(), regions[0, ])
  empty_cached <- helpers$formatSusieOutputs(list(), regions[0, ], extracted_results = list())
  for (name in table_names) {
    check(nrow(empty_tables[[name]]) == 0L &&
            identical(empty_tables[[name]]$feature_converged, logical()) &&
            identical(tail(names(empty_tables[[name]]), 1L), "feature_converged"),
          paste("Empty", name, "must append a logical convergence column"))
  }
  check(identical(empty_tables, empty_cached),
        "Empty cached extractions must retain the same convergence schema")
  empty_files <- helpers$writeSusieOutputs(list(), empty_tables,
                                          file.path(output_dir, "empty"))
  for (name in c("summary", "full")) {
    check(identical(arrow::read_parquet(empty_files[[name]])$feature_converged, logical()),
          paste("Empty serialized", name, "must retain a logical convergence column"))
  }

  # Keep the checkpointed formatter's shared constructor defaults unchanged.
  constructors <- list(helpers$InitEmptyCS(), helpers$InitEmptyInCSVariantDf(),
                       helpers$InitEmptyVariantDf())
  check(all(purrr::map_lgl(constructors, ~ !"feature_converged" %in% names(.x))),
        "Shared empty-table defaults must preserve the checkpointed output schema")
  if (length(failures)) stop(paste(failures, collapse = "\n"), call. = FALSE)

  message("SuSiE feature convergence tests passed: mixed outcomes, credible sets, cached output, serialization, old fits, and empty schemas")
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
  run_susie_feature_convergence_tests(repo_dir)
}
