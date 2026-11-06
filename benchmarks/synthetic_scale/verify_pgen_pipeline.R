# Compare complete results from the matched TSV and PGEN pipeline runs.
# Run this verifier separately from the measured workers.
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(readr))

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(1L, 2L)) {
  stop("Usage: Rscript verify_pgen_pipeline.R ROOT [LEGACY_PRODUCTION_DIR]",
       call. = FALSE)
}
root <- normalizePath(args[[1L]], mustWork = TRUE)
modes <- c("tabix_pipeline", "pgen_pipeline")
run_dirs <- stats::setNames(file.path(root, modes), modes)
stopifnot(all(dir.exists(run_dirs)))
checks <- list()
add_check <- function(name, passed, difference = NA_real_, detail = "") {
  checks[[length(checks) + 1L]] <<- tibble(
    check = name, passed = isTRUE(passed),
    max_abs_difference = difference, detail = detail
  )
}
numeric_difference <- function(x, y) {
  if (!is.numeric(x) || !is.numeric(y) || length(x) != length(y)) return(NA_real_)
  if (!identical(is.na(x), is.na(y)) || !identical(is.nan(x), is.nan(y))) return(Inf)
  finite <- is.finite(x) & is.finite(y)
  nonfinite <- !finite & !is.na(x)
  if (any(x[nonfinite] != y[nonfinite])) return(Inf)
  if (!any(finite)) return(0)
  max(abs(as.double(x[finite]) - as.double(y[finite])))
}
read_runs <- function(filename, reader) {
  paths <- file.path(run_dirs, filename)
  if (any(!file.exists(paths))) {
    stop("Missing benchmark output: ", paste(paths[!file.exists(paths)], collapse = ", "),
         call. = FALSE)
  }
  stats::setNames(map(paths, reader), modes)
}
compare_file <- function(filename, reader) {
  values <- read_runs(filename, reader)
  add_check(paste0(filename, "_exact"),
            identical(values$tabix_pipeline, values$pgen_pipeline))
  invisible(values)
}

fits <- compare_file("synthetic_susie.rds", readRDS)
tabix_fits <- fits$tabix_pipeline
pgen_fits <- fits$pgen_pipeline
stopifnot(is.list(tabix_fits), is.list(pgen_fits))
add_check("phenotype_names_and_order", identical(names(tabix_fits), names(pgen_fits)))
add_check("single_test_gene", length(tabix_fits) == 1L && length(pgen_fits) == 1L &&
            identical(names(tabix_fits), "synthetic_gene"))
walk(seq_len(max(length(tabix_fits), length(pgen_fits))), function(index) {
  if (index > length(tabix_fits) || index > length(pgen_fits)) {
    add_check(paste0("phenotype_", index, "_present"), FALSE)
    return(invisible(NULL))
  }
  old <- tabix_fits[[index]]
  new <- pgen_fits[[index]]
  label <- paste0("phenotype_", index, "_")
  add_check(paste0(label, "full_fit_exact"), identical(old, new))
  add_check(paste0(label, "PIPs_exact"), identical(old$pip, new$pip),
            numeric_difference(old$pip, new$pip))
  add_check(paste0(label, "credible_sets_exact"), identical(old$sets, new$sets))
  add_check(paste0(label, "variant_order_exact"), identical(old$variant_id, new$variant_id))
  add_check(paste0(label, "niter_exact"), identical(old$niter, new$niter),
            numeric_difference(old$niter, new$niter))
  add_check(paste0(label, "both_converged"), isTRUE(old$converged) && isTRUE(new$converged))
})
walk(c("synthetic.parquet", "synthetic.lbf_variable.parquet",
       "synthetic.full_susie.parquet"), function(filename) {
  compare_file(filename, ~ arrow::read_parquet(.x, as_data_frame = TRUE))
})
compare_file("synthetic.variant_position_summary.tsv",
             ~ read_lines(.x, progress = FALSE))
compare_file("synthetic.covariate_explained_variants.tsv",
             ~ read_lines(.x, progress = FALSE))
walk(c("rng_before.rds", "rng_after.rds"), ~ compare_file(.x, readRDS))
environments <- read_runs("environment.rds", readRDS)
environment_fields <- c("R", "BLAS", "susieR", "thread_variables")
add_check("same_R_BLAS_susieR_threads", identical(
  environments$tabix_pipeline[environment_fields],
  environments$pgen_pipeline[environment_fields]))

if (length(args) == 2L) {
  legacy_dir <- normalizePath(args[[2L]], mustWork = TRUE)
  legacy <- readRDS(file.path(legacy_dir, "synthetic_susie.rds"))
  old_fit_fields <- function(fit) {
    fit[setdiff(names(fit), "covariate_explained_variants")]
  }
  add_check("legacy_production_original_fit_fields_exact",
            identical(map(legacy, old_fit_fields), map(pgen_fits, old_fit_fields)),
            detail = legacy_dir)
}

timings <- imap_dfr(run_dirs, function(directory, run_mode) {
  usage <- read_tsv(file.path(directory, "resource_usage.tsv"), show_col_types = FALSE)
  stages <- read_tsv(file.path(directory, "stages.tsv"), show_col_types = FALSE)
  required_stages <- c("genotype_read", "fit", "entrypoint_total")
  add_check(paste0(run_mode, "_exit_zero"),
            nrow(usage) == 1L && identical(as.integer(usage$exit_code), 0L))
  add_check(paste0(run_mode, "_timings_complete"),
            all(required_stages %in% stages$stage))
  stage_seconds <- function(stage_name) {
    selected <- stages$elapsed_seconds[stages$stage == stage_name]
    if (length(selected) == 0L) return(NA_real_)
    sum(selected)
  }
  tibble(
    mode = run_mode,
    samples = usage$samples,
    variants = usage$variants,
    genotype_read_seconds = stage_seconds("genotype_read"),
    fit_seconds = stage_seconds("fit"),
    entrypoint_seconds = stage_seconds("entrypoint_total"),
    wall_seconds = usage$wall_seconds,
    peak_rss_gib = usage$peak_rss_gib
  )
})
result <- bind_rows(checks)
write_tsv(result, file.path(root, "pgen_checks.tsv"))
write_tsv(timings, file.path(root, "pgen_pipeline_times.tsv"))
print(result, n = Inf, width = Inf)
print(timings, width = Inf)
stopifnot(nrow(result) > 0L, all(result$passed))
cat("Matched TSV and PGEN pipeline outputs verified\n")
