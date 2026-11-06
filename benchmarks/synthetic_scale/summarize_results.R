suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(purrr))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
result_dir <- normalizePath(args[[1]], mustWork = TRUE)
modes <- c("baseline", "candidate", "factored")
read_stage <- function(stages, name) {
  values <- stages |> filter(stage == name) |> pull(elapsed_seconds)
  if (length(values) == 1L) values else NA_real_
}
metrics <- map_dfr(modes, function(mode) {
  directory <- file.path(result_dir, mode)
  usage <- read_tsv(file.path(directory, "resource_usage.tsv"), show_col_types = FALSE)
  stopifnot(nrow(usage) == 1L, usage$exit_code == 0L)
  stages <- read_tsv(file.path(directory, "stages.tsv"), show_col_types = FALSE)
  memory <- read_tsv(file.path(directory, "r_memory.tsv"), show_col_types = FALSE) |>
    filter(stage == "before_fit")
  comparison_file <- file.path(result_dir, paste0(mode, "_comparison"), "summary.json")
  comparison <- if (mode == "baseline") NULL else jsonlite::read_json(comparison_file, simplifyVector = TRUE)
  fits <- readRDS(file.path(directory, "synthetic_susie.rds"))
  fit <- fits[[1L]]
  tibble(mode = mode,
         exit_code = usage$exit_code,
         wall_seconds = usage$wall_seconds,
         entrypoint_seconds = read_stage(stages, "entrypoint_total"),
         entrypoint_without_fingerprints_seconds = read_stage(stages, "entrypoint_total") - read_stage(stages, "input_fingerprints"),
         genotype_residualization_seconds = read_stage(stages, "genotype_residualization"),
         phenotype_residualization_seconds = read_stage(stages, "phenotype_residualization"),
         fit_seconds = read_stage(stages, "fit"),
         output_processing_seconds = read_stage(stages, "output_processing"),
         fingerprint_seconds = read_stage(stages, "input_fingerprints"),
         peak_rss_GiB = usage$peak_rss_gib,
         r_used_before_fit_MiB = memory$r_used_MiB,
         variants = ncol(fit$alpha),
         credible_sets = length(fit$sets$cs),
         iterations = fit$niter,
         converged = fit$converged,
         exact_comparison_passed = if (mode == "baseline") NA else comparison$all_identical)
})
stopifnot(all(metrics$converged), all(metrics$variants == 20000L))
write_tsv(metrics, file.path(result_dir, "benchmark_results.tsv"))
baseline_fits <- readRDS(file.path(result_dir, "baseline", "synthetic_susie.rds"))
numerical_checks <- map_dfr(c("candidate", "factored"), function(mode) {
  fits <- readRDS(file.path(result_dir, mode, "synthetic_susie.rds"))
  map_dfr(names(baseline_fits), function(phenotype) {
    original <- baseline_fits[[phenotype]]
    changed <- fits[[phenotype]]
    tibble(mode = mode, phenotype = phenotype,
           complete_fit_identical = identical(original, changed),
           credible_set_membership_identical = identical(original$sets$cs, changed$sets$cs),
           max_pip_abs_difference = max(abs(original$pip - changed$pip)),
           max_alpha_abs_difference = max(abs(original$alpha - changed$alpha)),
           max_lbf_variable_abs_difference = max(abs(original$lbf_variable - changed$lbf_variable)),
           final_elbo_abs_difference = abs(tail(original$elbo, 1L) - tail(changed$elbo, 1L)),
           max_elbo_abs_difference = if (length(original$elbo) == length(changed$elbo)) {
             max(abs(original$elbo - changed$elbo))
           } else NA_real_)
  })
})
write_tsv(numerical_checks, file.path(result_dir, "numerical_checks.tsv"))
print(metrics, width = Inf)
print(numerical_checks, width = Inf)
stopifnot(isTRUE(metrics$exact_comparison_passed[metrics$mode == "candidate"]))
