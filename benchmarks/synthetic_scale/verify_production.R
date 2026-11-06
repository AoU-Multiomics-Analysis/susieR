suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[[1]], mustWork = TRUE)
read_fit <- function(mode) readRDS(file.path(root, mode, "synthetic_susie.rds"))[[1L]]
original <- read_fit("baseline")
accepted <- read_fit("factored")
production_with_diagnostics <- read_fit("production")
production <- production_with_diagnostics[
  setdiff(names(production_with_diagnostics), "covariate_explained_variants")
]
checks <- list()
add_check <- function(name, passed, difference = NA_real_) {
  checks[[length(checks) + 1L]] <<- tibble(check = name, passed = passed, max_abs_difference = difference)
}
usage <- read_tsv(file.path(root, "production", "resource_usage.tsv"), show_col_types = FALSE)
add_check("process_exit_zero", nrow(usage) == 1L && usage$exit_code == 0L)
diagnostics <- production_with_diagnostics$covariate_explained_variants
add_check("no_covariate_explained_variants_in_scale_fixture",
          is.data.frame(diagnostics) && nrow(diagnostics) == 0L &&
            identical(names(diagnostics), c("variant_id", "input_scale", "residual_sd",
                                           "residual_sd_ratio")))
add_check("same_fit_fields", identical(names(accepted), names(production)))
non_z <- setdiff(names(accepted), "z")
add_check("accepted_fit_except_z_exact", identical(accepted[non_z], production[non_z]))
add_check("PIPs_exact_to_accepted", identical(accepted$pip, production$pip),
          max(abs(accepted$pip - production$pip)))
add_check("credible_sets_match_original", identical(original$sets$cs, production$sets$cs))
add_check("all_variants_and_order", identical(original$variant_id, production$variant_id) &&
            ncol(production$alpha) == 20000L)
add_check("z_attributes", identical(attributes(accepted$z), attributes(production$z)))
z_difference <- max(abs(accepted$z - production$z))
add_check("z_agrees", is.finite(z_difference) && z_difference < 1e-10, z_difference)
for (filename in c("synthetic.parquet", "synthetic.lbf_variable.parquet", "synthetic.full_susie.parquet")) {
  old <- arrow::read_parquet(file.path(root, "factored", filename))
  new <- arrow::read_parquet(file.path(root, "production", filename))
  extra_fields <- setdiff(names(new), names(old))
  add_check(paste0(filename, "_diagnostic_fields"),
            length(extra_fields) == 0L ||
              identical(extra_fields, "has_covariate_explained_variant"))
  if ("has_covariate_explained_variant" %in% extra_fields) {
    add_check(paste0(filename, "_flags_false"),
              is.logical(new$has_covariate_explained_variant) &&
                all(!new$has_covariate_explained_variant))
  }
  new <- select(new, -any_of("has_covariate_explained_variant"))
  add_check(paste0(filename, "_original_fields"), identical(names(old), names(new)))
  add_check(paste0(filename, "_except_z"), identical(select(old, -any_of("z")), select(new, -any_of("z"))))
}
add_check("variant_summary_exact", identical(
  read_lines(file.path(root, "factored", "synthetic.variant_position_summary.tsv")),
  read_lines(file.path(root, "production", "synthetic.variant_position_summary.tsv"))))
for (filename in c("rng_before.rds", "rng_after.rds")) {
  add_check(filename, identical(readRDS(file.path(root, "factored", filename)),
                                readRDS(file.path(root, "production", filename))))
}
result <- bind_rows(checks)
write_tsv(result, file.path(root, "production_checks.tsv"))
stages <- read_tsv(file.path(root, "production", "stages.tsv"), show_col_types = FALSE)
metrics <- tibble(
  samples = 10000L, variants = ncol(production$alpha),
  entrypoint_seconds = stages$elapsed_seconds[stages$stage == "entrypoint_total"],
  wall_seconds = usage$wall_seconds, peak_rss_GiB = usage$peak_rss_gib,
  iterations = production$niter, converged = production$converged,
  credible_sets = length(production$sets$cs),
  max_PIP_difference_from_original = max(abs(production$pip - original$pip)),
  max_PIP_difference_from_accepted = max(abs(production$pip - accepted$pip)),
  max_z_difference_from_accepted = z_difference)
write_tsv(metrics, file.path(root, "production_results.tsv"))
print(result, n = Inf)
print(metrics, width = Inf)
stopifnot(all(result$passed), production$converged)
cat("Actual production run verified\n")
