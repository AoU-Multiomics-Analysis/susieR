.libPaths(c(Sys.getenv("BENCH_LIBRARY"), .libPaths()))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 6L)
input_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
mode <- match.arg(args[[3]], c("baseline", "candidate", "factored", "factored_fast"))
repo_dir <- normalizePath(args[[4]], mustWork = TRUE)
n_samples <- as.integer(args[[5]])
n_variants <- as.integer(args[[6]])
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)
stopifnot(requireNamespace("eQTLUtils", quietly = TRUE))

bench_times <- list()
bench_stages <- tibble(stage = character(), elapsed_seconds = double())
bench_memory <- tibble(stage = character(), r_used_MiB = double(), r_vector_used_MiB = double())
bench_begin <- function(stage) {
  bench_times[[stage]] <<- proc.time()[["elapsed"]]
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " START ", stage)
}
bench_end <- function(stage) {
  elapsed <- proc.time()[["elapsed"]] - bench_times[[stage]]
  bench_stages <<- bind_rows(bench_stages, tibble(stage = stage, elapsed_seconds = elapsed))
  write_tsv(bench_stages, file.path(output_dir, "stages.tsv"))
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " END ", stage, ": ", round(elapsed, 3), " s")
}
bench_gc <- function(stage) {
  stats <- gc()
  bench_memory <<- bind_rows(bench_memory, tibble(stage = stage,
                                                r_used_MiB = sum(stats[, 2L]),
                                                r_vector_used_MiB = stats["Vcells", 2L]))
  write_tsv(bench_memory, file.path(output_dir, "r_memory.tsv"))
  invisible(stats)
}
bench_hash_matrix <- function(mat) {
  column_hashes <- vapply(seq_len(ncol(mat)), function(j) {
    digest::digest(mat[, j], algo = "sha256")
  }, character(1))
  digest::digest(list(type = typeof(mat), dimensions = dim(mat), names = dimnames(mat),
                      columns = column_hashes), algo = "sha256")
}
bench_before_fit <- function(gt_hat, expression_vector) {
  stopifnot(identical(dim(gt_hat), c(n_samples, n_variants)), all(is.finite(expression_vector)))
  bench_begin("input_fingerprints")
  fingerprints <- tibble(object = c("residualized_genotypes", "residualized_phenotype"),
                          sha256 = c(bench_hash_matrix(gt_hat), digest::digest(expression_vector, algo = "sha256")))
  write_tsv(fingerprints, file.path(output_dir, "input_fingerprints.tsv"))
  bench_end("input_fingerprints")
  saveRDS(.Random.seed, file.path(output_dir, "rng_before.rds"))
  cache_path <- Sys.getenv("BENCH_FIT_CACHE", unset = "")
  if (nzchar(cache_path)) {
    message("Save residualized fit inputs to temporary cache")
    saveRDS(list(gt_hat = gt_hat, expression_vector = expression_vector,
                 fit_seed = .Random.seed, environment = environment_info),
            cache_path, compress = FALSE)
    if (identical(Sys.getenv("BENCH_PREPARE_ONLY"), "true")) {
      message("Fit-input preparation completed")
      quit(save = "no", status = 0L)
    }
  }
  bench_begin("fit")
}
bench_after_fit <- function() {
  bench_end("fit")
  saveRDS(.Random.seed, file.path(output_dir, "rng_after.rds"))
}
replace_once <- function(text, old, new) {
  matches <- gregexpr(old, text, fixed = TRUE)[[1L]]
  stopifnot(length(matches) == 1L, matches[[1L]] > 0L)
  sub(old, new, text, fixed = TRUE)
}

helper_dir <- file.path(output_dir, "helpers")
dir.create(helper_dir, showWarnings = FALSE)
helper_paths <- list.files(file.path(repo_dir, "R", "utils"), pattern = "[.]R$", full.names = TRUE)
stopifnot(all(file.copy(helper_paths, helper_dir, overwrite = TRUE)))
helper_path <- file.path(helper_dir, "SusieFunctions.R")
helper <- paste(readLines(helper_path), collapse = "\n")
if (mode != "baseline") {
  helper <- replace_once(helper,
    "gt_hat <- hat %*% t(gt_std)  # result: variants x samples\n  rm(gt_std)",
    "gt_hat <- hat %*% t(gt_std)  # unchanged arithmetic\n  rm(gt_std, hat)")
  begin <- regexpr("extractPipsFromVariantDf <- function(variant_df){", helper, fixed = TRUE)[[1L]]
  end <- regexpr("# Build connected components across credible sets", helper, fixed = TRUE)[[1L]]
  stopifnot(begin > 0L, end > begin)
  replacement <- paste(c(
    "extractPipsFromVariantDf <- function(variant_df) {",
    "  rows <- variant_df %>% dplyr::filter(!is.na(cs_index))",
    "  keys <- rows %>% dplyr::select(phenotype_id, variant_id, cs_id, cs_index)",
    "  if (nrow(rows) == 0L) return(dplyr::mutate(keys, pip = rows$alpha1))",
    "  alpha <- rows %>% dplyr::select(alpha1:alpha10) %>% as.matrix()",
    "  component <- match(rows$cs_index, paste0('L', seq_len(10L)))",
    "  keep <- which(!is.na(component))",
    "  keys %>% dplyr::slice(keep) %>% dplyr::mutate(pip = alpha[cbind(keep, component[keep])])",
    "}", "", ""), collapse = "\n")
  helper <- paste0(substr(helper, 1L, begin - 1L), replacement, substr(helper, end, nchar(helper)))
}
if (mode %in% c("factored", "factored_fast")) {
  helper <- replace_once(helper,
    "hat = diag(nrow(covariates_matrix)) - covariates_matrix %*% solve(crossprod(covariates_matrix)) %*% t(covariates_matrix)\n  expression_vector = hat %*% gene_vector$phenotype_value_std",
    "projection <- covariates_matrix %*% solve(crossprod(covariates_matrix))\n  expression_vector <- gene_vector$phenotype_value_std - projection %*% crossprod(covariates_matrix, gene_vector$phenotype_value_std)")
  helper <- replace_once(helper,
    "gt_hat <- hat %*% t(gt_std)  # unchanged arithmetic\n  rm(gt_std, hat)",
    "gt_transpose <- t(gt_std)\n  gt_hat <- gt_transpose - projection %*% crossprod(covariates_matrix, gt_transpose)\n  rm(gt_transpose, gt_std, projection)")
}
helper <- replace_once(helper, "message('Residualizing gene expression')",
                       "bench_begin('phenotype_residualization')\n  message('Residualizing gene expression')")
helper <- replace_once(helper, "names(expression_vector) = gene_vector$genotype_id",
                       "names(expression_vector) = gene_vector$genotype_id\n  bench_end('phenotype_residualization')")
helper <- replace_once(helper, "message('Computing gt_hat')",
                       "bench_begin('genotype_residualization')\n  message('Computing gt_hat')")
if (mode == "baseline") {
  helper <- replace_once(helper,
    "gt_hat <- hat %*% t(gt_std)  # result: variants x samples\n  rm(gt_std)\n  gc()",
    "gt_hat <- hat %*% t(gt_std)  # result: variants x samples\n  bench_end('genotype_residualization')\n  rm(gt_std)\n  bench_gc('before_fit')")
} else if (mode == "candidate") {
  helper <- replace_once(helper,
    "gt_hat <- hat %*% t(gt_std)  # unchanged arithmetic\n  rm(gt_std, hat)\n  gc()",
    "gt_hat <- hat %*% t(gt_std)  # unchanged arithmetic\n  bench_end('genotype_residualization')\n  rm(gt_std, hat)\n  bench_gc('before_fit')")
} else {
  helper <- replace_once(helper,
    "gt_hat <- gt_transpose - projection %*% crossprod(covariates_matrix, gt_transpose)\n  rm(gt_transpose, gt_std, projection)\n  gc()",
    "gt_hat <- gt_transpose - projection %*% crossprod(covariates_matrix, gt_transpose)\n  bench_end('genotype_residualization')\n  rm(gt_transpose, gt_std, projection)\n  bench_gc('before_fit')")
}
helper <- replace_once(helper,
  "# Fit the SuSiE model on residualized genotypes and expression.\n  fitted <- susieR::susie",
  "# Fit the SuSiE model on residualized genotypes and expression.\n  bench_before_fit(gt_hat, expression_vector)\n  fitted <- susieR::susie")
helper <- replace_once(helper,
  "min_abs_corr = 0.5)\n  fitted$variant_id = variant_names",
  "min_abs_corr = 0.5)\n  bench_after_fit()\n  fitted$variant_id = variant_names")
if (mode == "candidate") {
  helper <- replace_once(helper,
    "gt_hat <- hat %*% t(gt_std)  # unchanged arithmetic",
    "gt_hat <- tcrossprod(hat, gt_std)  # avoid the explicit transpose")
}
if (mode == "factored_fast") {
  for (prototype in c("fast_zscore.R", "fast_covariates.R")) {
    snapshot <- file.path(output_dir, prototype)
    stopifnot(file.copy(file.path(repo_dir, "benchmarks", "synthetic_scale", prototype),
                       snapshot, overwrite = TRUE))
    source(snapshot)
  }
  helper <- replace_once(helper,
    "verbose = TRUE,\n                          compute_univariate_zscore = TRUE,",
    "verbose = TRUE,\n                          compute_univariate_zscore = FALSE,")
  helper <- replace_once(helper,
    "bench_after_fit()\n  fitted$variant_id = variant_names",
    "fitted$z <- fast_univariate_zscore(gt_hat, expression_vector)\n  bench_after_fit()\n  fitted$variant_id = variant_names")
  import_path <- file.path(helper_dir, "ImportFunctions.R")
  import_helper <- paste(readLines(import_path), collapse = "\n")
  import_helper <- replace_once(import_helper,
    "pc_matrix = read.table(covariates_path, check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)",
    "pc_matrix = fast_read_covariates(covariates_path)")
  write_lines(import_helper, import_path)
}
write_lines(helper, helper_path)

entrypoint <- paste(readLines(file.path(repo_dir, "R", "scripts", "susie.R")), collapse = "\n")
entrypoint <- replace_once(entrypoint,
  "opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))",
  "opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list), args = bench_cli_args)")
entrypoint <- replace_once(entrypoint, "data_list <- LoadData(opt)",
                           "bench_begin('load_data')\ndata_list <- LoadData(opt)\nbench_end('load_data')")
entrypoint <- replace_once(entrypoint, "results = purrr::map(selected_phenotypes,",
                           "bench_begin('fine_mapping_total')\nresults = purrr::map(selected_phenotypes,")
entrypoint <- replace_once(entrypoint, "saveRDS(results,file = paste0(opt$out_prefix,'_susie.rds') )",
                           "bench_end('fine_mapping_total')\nsaveRDS(results,file = paste0(opt$out_prefix,'_susie.rds') )\nbench_begin('output_processing')")
entrypoint <- replace_once(entrypoint, "quit(save = \"no\", status = 0)",
                           "bench_end('output_processing')\n  bench_end('entrypoint_total')\n  quit(save = \"no\", status = 0)")
if (mode != "baseline") {
  entrypoint <- replace_once(entrypoint,
    "res = purrr::map(results, extractResults) %>%\n    purrr::transpose()",
    "res = purrr::map(results, extractResults) %>%\n    purrr::transpose()\nrm(results, reusable_genotype)\ninvisible(gc())")
  entrypoint <- replace_once(entrypoint,
    "in_cs_variant_gene_df <- in_cs_variant_df %>% \n  dplyr::left_join(phenotype_meta %>% dplyr::select(phenotype_id, gene_id, group_id), by = c(\"molecular_trait_id\" = \"phenotype_id\")) %>% \n  dplyr::group_by(gene_id) %>% \n  dplyr::mutate(uniq_phenotypes_count = length(base::unique(molecular_trait_id))) %>% \n  dplyr::ungroup()",
    "# The unused connected-component table is omitted in this benchmark candidate.")
}
entrypoint <- paste0(entrypoint, "\nbench_end('output_processing')\n")
entrypoint_path <- file.path(output_dir, "entrypoint.R")
write_lines(entrypoint, entrypoint_path)
Sys.setenv(SUSIER_FUNCTIONS_PATH = helper_dir)
bench_cli_args <- c("--genotype_matrix", file.path(input_dir, "dosages.tsv.bgz"),
                    "--sample_meta", file.path(input_dir, "samples.tsv"),
                    "--phenotype_list", file.path(input_dir, "permutations.tsv"),
                    "--expression_matrix", file.path(input_dir, "phenotypes.bed"),
                    "--covariates", file.path(input_dir, "covariates.tsv"),
                    "--out_prefix", file.path(output_dir, "synthetic"),
                    "--cisdistance", "1000000", "--MAF", "0")
versions <- c("susieR", "eQTLUtils", "dplyr", "tidyr", "readr", "arrow", "Rsamtools")
environment_info <- list(R = R.version.string, platform = R.version$platform,
                          package_versions = setNames(vapply(versions, function(p) as.character(packageVersion(p)), character(1)), versions),
                          BLAS = extSoftVersion()[["BLAS"]],
                          thread_variables = Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "ARROW_NUM_THREADS")),
                          eQTLUtils_revision = "26242562a4e244334fd9691d03bc1ef4d2d6c1d9",
                          dimensions = c(samples = n_samples, variants = n_variants),
                          seed = 20261106L)
saveRDS(environment_info, file.path(output_dir, "environment.rds"))
write_lines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
set.seed(20261106)
bench_begin("entrypoint_total")
source(entrypoint_path, local = .GlobalEnv)
bench_end("entrypoint_total")
message("Benchmark worker completed: ", mode)
