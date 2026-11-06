suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(purrr))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[[1]], mustWork = TRUE)
read_fit <- function(mode) readRDS(file.path(root, mode, "synthetic_susie.rds"))[[1L]]
reference <- read_fit("factored")
combined <- read_fit("factored_fast")
check_rows <- list()
add_check <- function(check, passed, max_difference = NA_real_) {
  check_rows[[length(check_rows) + 1L]] <<- tibble(check = check, passed = passed,
                                                max_abs_difference = max_difference)
}
non_z <- setdiff(names(reference), "z")
add_check("full_fit_fields", identical(names(reference), names(combined)))
add_check("full_fit_except_z", identical(reference[non_z], combined[non_z]))
add_check("z_attributes", identical(attributes(reference$z), attributes(combined$z)))
z_difference <- max(abs(reference$z - combined$z))
add_check("z_numerical_agreement", is.finite(z_difference) && z_difference < 1e-10,
          z_difference)
add_check("PIPs_exact", identical(reference$pip, combined$pip),
          max(abs(reference$pip - combined$pip)))
add_check("credible_set_membership", identical(reference$sets$cs, combined$sets$cs))
for (filename in c("synthetic.parquet", "synthetic.lbf_variable.parquet",
                   "synthetic.full_susie.parquet")) {
  old <- arrow::read_parquet(file.path(root, "factored", filename))
  new <- arrow::read_parquet(file.path(root, "factored_fast", filename))
  add_check(paste0(filename, "_fields"), identical(names(old), names(new)))
  add_check(paste0(filename, "_except_z"), identical(select(old, -any_of("z")),
                                                   select(new, -any_of("z"))))
}
for (filename in c("synthetic.variant_position_summary.tsv", "input_fingerprints.tsv")) {
  add_check(filename, identical(read_lines(file.path(root, "factored", filename)),
                                 read_lines(file.path(root, "factored_fast", filename))))
}
for (filename in c("rng_before.rds", "rng_after.rds", "environment.rds")) {
  add_check(filename, identical(readRDS(file.path(root, "factored", filename)),
                                 readRDS(file.path(root, "factored_fast", filename))))
}

metrics <- map_dfr(c("baseline", "factored", "factored_fast"), function(mode) {
  usage <- read_tsv(file.path(root, mode, "resource_usage.tsv"), show_col_types = FALSE)
  stages <- read_tsv(file.path(root, mode, "stages.tsv"), show_col_types = FALSE)
  stopifnot(nrow(usage) == 1L, usage$exit_code == 0L)
  stage_time <- function(stage_name) {
    result <- stages |> filter(stage == stage_name) |> pull(elapsed_seconds)
    stopifnot(length(result) == 1L)
    result
  }
  fit <- read_fit(mode)
  stopifnot(isTRUE(fit$converged), ncol(fit$alpha) == 20000L)
  tibble(mode = mode, wall_seconds = usage$wall_seconds,
         run_without_fingerprints_seconds = stage_time("entrypoint_total") - stage_time("input_fingerprints"),
         load_data_seconds = stage_time("load_data"),
         fit_and_z_seconds = stage_time("fit"), peak_rss_GiB = usage$peak_rss_gib,
         iterations = fit$niter, credible_sets = length(fit$sets$cs))
})
fit_metrics <- map_dfr(c("fit_quiet", "fit_fast_z", "fit_fast_z_threads4"), function(mode) {
  timings <- read_tsv(file.path(root, mode, "fit_timings.tsv"), show_col_types = FALSE)
  usage <- read_tsv(file.path(root, mode, "resource_usage.tsv"), show_col_types = FALSE)
  stopifnot(usage$exit_code == 0L, timings$converged)
  fit <- read_fit(mode)
  tibble(mode = mode, threads = usage$threads_requested,
         fit_and_z_seconds = timings$fit_and_z_seconds,
         PIPs_exact = identical(reference$pip, fit$pip),
         credible_set_membership_exact = identical(reference$sets$cs, fit$sets$cs),
         model_except_z_exact = identical(reference[non_z], fit[non_z]),
         max_z_abs_difference = max(abs(reference$z - fit$z)))
})
checks <- bind_rows(check_rows)
write_tsv(checks, file.path(root, "followup_checks.tsv"))
write_tsv(metrics, file.path(root, "followup_results.tsv"))
write_tsv(fit_metrics, file.path(root, "followup_fit_results.tsv"))
print(metrics, width = Inf)
print(fit_metrics, width = Inf)
stopifnot(all(checks$passed), all(fit_metrics$PIPs_exact),
          all(fit_metrics$credible_set_membership_exact))
cat("Follow-up validation passed\n")
