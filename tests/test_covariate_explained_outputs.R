# Exercise the real entry point's export expressions with known fitted results.
# Removing the TSV export or the all-members CS flag must fail these checks.
run_covariate_explained_output_tests <- function(repo_dir = ".") {
  functions_dir <- Sys.getenv("SUSIER_FUNCTIONS_PATH",
                             unset = file.path(repo_dir, "R", "utils"))
  entrypoint_paths <- c(file.path(repo_dir, "R", "scripts", "susie.R"),
                       file.path(dirname(functions_dir), "scripts", "susie.R"))
  entrypoint <- entrypoint_paths[file.exists(entrypoint_paths)][[1L]]
  test_dir <- tempfile("susier-diagnostic-outputs-")
  dir.create(test_dir)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  set.seed(479L)
  variant_ids <- paste0("chr1_", 1000L + seq_len(12L), "_A_G")
  x <- matrix(rnorm(100L * 12L), nrow = 100L,
              dimnames = list(NULL, variant_ids))
  y <- cbind(2 * x[, 1L] + x[, 3L] + rnorm(100L))
  fit <- susieR::susie(x, y, L = 10, compute_univariate_zscore = TRUE)
  fit$variant_id <- variant_ids
  # Fix set membership for output tests; retain all fitted numeric values.
  fit$sets$cs <- list(L1 = c(1L, 2L), L2 = c(3L, 4L))
  fit$sets$cs_index <- c(1L, 2L)
  fit$sets$purity <- data.frame(
    min.abs.corr = c(0.9, 0.8), mean.abs.corr = c(0.95, 0.9),
    median.abs.corr = c(0.95, 0.9), row.names = c("L1", "L2")
  )
  diagnostic <- tibble::tibble(
    variant_id = variant_ids[c(2L, 9L)], input_scale = c(1, 2),
    residual_sd = c(1e-10, 4e-10), residual_sd_ratio = c(1e-10, 2e-10)
  )
  annotated <- fit
  annotated$covariate_explained_variants <- diagnostic
  empty_diagnostic <- diagnostic[0, ]
  no_flags <- fit
  no_flags$covariate_explained_variants <- empty_diagnostic
  no_sets <- annotated
  no_sets$sets$cs <- NULL
  no_sets$sets$cs_index <- integer()
  no_sets$sets$purity <- fit$sets$purity[0, , drop = FALSE]

  # Run only the entry point's output stage. Model fitting/input loading have
  # separate tests. A child process exercises the shared formatter and writer calls.
  driver <- file.path(test_dir, "export_driver.R")
  writeLines(c(
    "args <- commandArgs(trailingOnly = TRUE)",
    "bench_library <- Sys.getenv('BENCH_LIBRARY', unset = '')",
    "if (nzchar(bench_library)) .libPaths(c(bench_library, .libPaths()))",
    "suppressPackageStartupMessages(library(dplyr))",
    "suppressPackageStartupMessages(library(tidyr))",
    "suppressPackageStartupMessages(library(stringr))",
    "suppressPackageStartupMessages(library(purrr))",
    "source(file.path(args[[1L]], 'InitFunctions.R'))",
    "source(file.path(args[[1L]], 'SusieFunctions.R'))",
    "source(file.path(args[[1L]], 'SusieOutputFunctions.R'))",
    "results <- readRDS(args[[3L]])",
    "opt <- list(out_prefix = args[[4L]])",
    "region_df <- tibble(phenotype_id = names(results), region = 'chr1_1_2000')",
    "phenotype_meta <- tibble(phenotype_id = names(results), gene_id = names(results), group_id = names(results))",
    "phenotype_list <- phenotype_meta",
    "selected_group_ids <- names(results)",
    "selected_phenotypes <- as.list(names(results))",
    "expressions <- parse(file = args[[2L]])",
    "for (expression in expressions) {",
    "  if (is.call(expression) && identical(expression[[1L]], as.name('<-')) &&",
    "      grepl('^empty_', as.character(expression[[2L]]))) eval(expression)",
    "}",
    "first_export <- which(vapply(expressions, function(expression) {",
    "  is.call(expression) && identical(expression[[1L]], as.name('<-')) &&",
    "    identical(expression[[2L]], as.name('tables'))",
    "}, logical(1)))[[1L]]",
    "eval(as.expression(expressions[seq.int(first_export, length(expressions))]))",
    "saveRDS(tables$cs_df, paste0(opt$out_prefix, '.cs_summary.rds'))"
  ), driver)
  run_export <- function(results, label) {
    fixture_path <- file.path(test_dir, paste0(label, ".rds"))
    saveRDS(results, fixture_path)
    prefix <- file.path(test_dir, label)
    log_path <- paste0(prefix, ".log")
    status <- system2(file.path(R.home("bin"), "Rscript"),
      args = shQuote(c(driver, normalizePath(functions_dir),
                       normalizePath(entrypoint), fixture_path, prefix)),
      stdout = log_path, stderr = log_path)
    if (!identical(status, 0L)) {
      stop(label, " export failed:\n", paste(readLines(log_path), collapse = "\n"),
           call. = FALSE)
    }
    diagnostic_path <- paste0(prefix, ".covariate_explained_variants.tsv")
    if (!file.exists(diagnostic_path)) {
      stop("The entry point must always write the covariate-explained variant TSV.",
           call. = FALSE)
    }
    list(
      diagnostic = readr::read_tsv(diagnostic_path, show_col_types = FALSE,
        col_types = "ccddd"),
      cs = arrow::read_parquet(paste0(prefix, ".parquet")),
      full = arrow::read_parquet(paste0(prefix, ".full_susie.parquet")),
      summary = if (file.exists(paste0(prefix, ".cs_summary.rds"))) {
        readRDS(paste0(prefix, ".cs_summary.rds"))
      } else NULL
    )
  }

  flagged <- run_export(list(gene1 = annotated, gene2 = no_flags), "flagged")
  expected_diagnostic <- diagnostic |>
    dplyr::transmute(molecular_trait_id = "gene1", variant = variant_id,
                     input_scale, residual_sd, residual_sd_ratio)
  stopifnot(identical(as.data.frame(flagged$diagnostic),
                      as.data.frame(expected_diagnostic)))
  stopifnot(is.logical(flagged$cs$has_covariate_explained_variant),
            is.logical(flagged$summary$has_covariate_explained_variant))
  expected_flag <- flagged$cs$molecular_trait_id == "gene1" &
    flagged$cs$cs_index == "L1"
  stopifnot(identical(flagged$cs$has_covariate_explained_variant, expected_flag),
            all(flagged$cs$has_covariate_explained_variant[
              flagged$cs$molecular_trait_id == "gene1" & flagged$cs$cs_index == "L1"]),
            identical(flagged$summary$has_covariate_explained_variant,
                      flagged$summary$molecular_trait_id == "gene1" &
                        flagged$summary$cs_index == "L1"))

  old <- run_export(list(gene1 = fit, gene2 = fit), "old_fit")
  stopifnot(nrow(old$diagnostic) == 0L,
            identical(names(old$diagnostic), names(expected_diagnostic)),
            !any(old$cs$has_covariate_explained_variant),
            !any(old$summary$has_covariate_explained_variant))
  numeric_columns <- names(old$cs)[vapply(old$cs, is.numeric, logical(1))]
  stopifnot(identical(as.data.frame(old$cs[numeric_columns]),
                      as.data.frame(flagged$cs[numeric_columns])),
            identical(as.data.frame(old$full), as.data.frame(flagged$full)))

  zero_flags <- run_export(list(gene1 = no_flags), "zero_flags")
  stopifnot(nrow(zero_flags$diagnostic) == 0L,
            !any(zero_flags$cs$has_covariate_explained_variant))

  no_credible_sets <- run_export(list(gene1 = no_sets), "no_credible_sets")
  stopifnot(nrow(no_credible_sets$diagnostic) == nrow(diagnostic),
            nrow(no_credible_sets$cs) == 0L,
            is.logical(no_credible_sets$cs$has_covariate_explained_variant))
  no_sets$covariate_explained_variants <- empty_diagnostic
  empty <- run_export(list(gene1 = no_sets), "empty_diagnostic")
  stopifnot(nrow(empty$diagnostic) == 0L,
            identical(names(empty$diagnostic), names(expected_diagnostic)),
            nrow(empty$cs) == 0L,
            is.logical(empty$cs$has_covariate_explained_variant))
  message("Covariate-explained output tests passed: TSV, CS flags, old fits, and no-set outputs")
  invisible(TRUE)
}
