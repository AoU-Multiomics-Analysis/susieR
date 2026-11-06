# Isolated model benchmark. The same cached X, y, and RNG state feed every run.
.libPaths(c(Sys.getenv("BENCH_LIBRARY"), .libPaths()))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(susieR))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 6L)
cache_file <- args[[1]]
output_dir <- args[[2]]
mode <- match.arg(args[[3]], c("fit_baseline", "fit_quiet", "fit_fast_z"))
repo_dir <- args[[4]]
n_samples <- as.integer(args[[5]])
n_variants <- as.integer(args[[6]])
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
message("Load the shared residualized input cache")
load_time <- system.time(cache <- readRDS(cache_file))[["elapsed"]]
gt_hat <- cache$gt_hat
expression_vector <- cache$expression_vector
fit_seed <- cache$fit_seed
cache_environment <- cache$environment
rm(cache)
invisible(gc())
stopifnot(identical(dim(gt_hat), c(n_samples, n_variants)))
if (mode == "fit_fast_z") {
  z_snapshot <- file.path(output_dir, "fast_zscore_snapshot.R")
  stopifnot(file.copy(file.path(repo_dir, "benchmarks", "synthetic_scale", "fast_zscore.R"),
                     z_snapshot, overwrite = TRUE))
  source(z_snapshot)
}
saveRDS(list(BLAS = extSoftVersion()[["BLAS"]],
             R = R.version.string,
             susieR = as.character(packageVersion("susieR")),
             threads = Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")),
             cache_environment = cache_environment),
        file.path(output_dir, "environment.rds"))
assign(".Random.seed", fit_seed, envir = .GlobalEnv)
saveRDS(.Random.seed, file.path(output_dir, "rng_before.rds"))
message("Fit the same 10-effect SuSiE model: ", mode)
fit_time <- system.time({
  fitted <- susieR::susie(gt_hat, expression_vector,
                         L = 10,
                         estimate_residual_variance = TRUE,
                         estimate_prior_variance = TRUE,
                         scaled_prior_variance = 0.1,
                         verbose = mode == "fit_baseline",
                         compute_univariate_zscore = mode != "fit_fast_z",
                         min_abs_corr = 0.5)
})[["elapsed"]]
z_time <- 0
if (mode == "fit_fast_z") {
  message("Calculate marginal z-scores with the vectorized method")
  z_time <- system.time({
    fitted$z <- fast_univariate_zscore(gt_hat, expression_vector)
  })[["elapsed"]]
}
fitted$variant_id <- colnames(gt_hat)
saveRDS(list(synthetic_gene = fitted), file.path(output_dir, "synthetic_susie.rds"))
saveRDS(.Random.seed, file.path(output_dir, "rng_after.rds"))
write_tsv(tibble(mode = mode, cache_load_seconds = load_time,
                 susie_seconds = fit_time, added_z_seconds = z_time,
                 fit_and_z_seconds = fit_time + z_time,
                 variants = ncol(fitted$alpha), credible_sets = length(fitted$sets$cs),
                 iterations = fitted$niter, converged = fitted$converged),
          file.path(output_dir, "fit_timings.tsv"))
message("Isolated fit completed: ", fit_time + z_time, " s")
