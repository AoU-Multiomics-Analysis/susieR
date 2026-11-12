run_finemap_path_tests <- function(repo_dir = ".") {
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(tidyr))
  suppressPackageStartupMessages(library(purrr))
  functions_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH", unset = file.path(repo_dir, "R", "utils"))
  for (filename in c("SusiePerformance.R", "ImportFunctions.R", "MemoryHelpers.R", "SusieFunctions.R")) {
    if (file.exists(file.path(functions_path, filename))) source(file.path(functions_path, filename), local = environment())
  }
  assert_flag_table <- function(table, variant_ids = character()) {
    stopifnot("Covariate-explained diagnostics must be a typed tibble" =
                tibble::is_tibble(table),
              identical(names(table), c("variant_id", "input_scale", "residual_sd", "residual_sd_ratio")),
              identical(vapply(table, typeof, character(1)),
                        c(variant_id = "character", input_scale = "double",
                          residual_sd = "double", residual_sd_ratio = "double")),
              identical(table$variant_id, variant_ids),
              all(table$input_scale > 0), all(table$residual_sd_ratio <= 1e-8),
              identical(table$residual_sd_ratio, table$residual_sd / table$input_scale))
  }

  # A covariate-explained predictor has only numerical residuals. SuSiE can
  # amplify a change in those residuals when it standardizes the predictor.
  # Keep the complete old matrix multiplication as an independent reference.
  set.seed(11L)
  guard_n <- 180L
  guard_ids <- paste0("g", seq_len(guard_n))
  sex <- rbinom(guard_n, 1L, 0.5)
  guard_design <- cbind(sex = sex, pc = rnorm(guard_n), intercept = 1)
  rownames(guard_design) <- guard_ids
  guard_raw_x <- cbind(sex = sex, matrix(rbinom(guard_n * 15L, 2L, 0.3), nrow = guard_n))
  guard_x <- guard_raw_x
  guard_x <- sweep(guard_x, 2L, colMeans(guard_x), "-")
  rownames(guard_x) <- guard_ids
  colnames(guard_x) <- paste0("variant", seq_len(ncol(guard_x)))
  guard_y <- rnorm(guard_n)
  guard_y <- qnorm((rank(guard_y) - 0.5) / guard_n)
  names(guard_y) <- guard_ids
  guard_adjustment <- createCovariateAdjustment(guard_design)
  guard_hat <- diag(guard_n) - guard_design %*% solve(crossprod(guard_design)) %*% t(guard_design)
  legacy_x <- guard_hat %*% guard_x
  factored_x <- guard_x - guard_adjustment$projection %*% crossprod(guard_design, guard_x)
  guard_rng <- .Random.seed
  guarded_x <- residualizeWithCovariates(guard_adjustment, guard_x,
                                        variant_ids = colnames(guard_x))
  stopifnot("Explained predictors must use the complete legacy multiplication" =
              identical(guarded_x, legacy_x),
            "The explained fixture must distinguish the two calculations" =
              !identical(factored_x, legacy_x),
            identical(guard_rng, .Random.seed))
  assert_flag_table(guard_adjustment$state$covariate_explained_variants, colnames(guard_x)[1L])
  guard_diagnostics <- guard_adjustment$state$covariate_explained_variants
  invisible(residualizeWithCovariates(guard_adjustment, guard_y))
  stopifnot("Response calls must retain genotype diagnostics" =
              identical(guard_adjustment$state$covariate_explained_variants, guard_diagnostics))

  # A zero column must not send unrelated, well-behaved predictors to the
  # dense path. Their numerical results must remain unchanged.
  ordinary_x <- cbind(zero = 0, guard_x[, -1L, drop = FALSE])
  ordinary_adjustment <- createCovariateAdjustment(guard_design)
  ordinary_factored <- ordinary_x - guard_adjustment$projection %*% crossprod(guard_design, ordinary_x)
  stopifnot("Ordinary predictors and zero columns must retain the fast result" =
              identical(residualizeWithCovariates(ordinary_adjustment, ordinary_x,
                                                  variant_ids = colnames(ordinary_x)), ordinary_factored),
            "The ordinary fixture must distinguish the two calculations" =
              !identical(ordinary_factored, guard_hat %*% ordinary_x))
  assert_flag_table(ordinary_adjustment$state$covariate_explained_variants)
  nonzero_constant <- cbind(constant = 1, ordinary_x)
  stopifnot("An explained nonzero constant must retain its legacy residuals" =
              identical(residualizeWithCovariates(guard_adjustment, nonzero_constant),
                        guard_hat %*% nonzero_constant))

  # The same guard must apply to a response vector whose covariate contains
  # the rank-normalized response plus a small perturbation.
  near_design <- cbind(pc = guard_y + 1e-10 * rnorm(guard_n), intercept = 1)
  rownames(near_design) <- guard_ids
  near_adjustment <- createCovariateAdjustment(near_design)
  near_hat <- diag(guard_n) - near_design %*% solve(crossprod(near_design)) %*% t(near_design)
  legacy_y <- near_hat %*% guard_y
  guarded_y <- residualizeWithCovariates(near_adjustment, guard_y)
  stopifnot("A nearly explained response must retain the legacy residuals" =
              identical(guarded_y, legacy_y),
            identical(dim(guarded_y), c(guard_n, 1L)),
            identical(rownames(guarded_y), guard_ids))

  legacy_fit <- function(x, y) {
    susieR::susie(x, y, L = 10,
                  estimate_residual_variance = TRUE,
                  estimate_prior_variance = TRUE,
                  scaled_prior_variance = 0.1, verbose = FALSE,
                  compute_univariate_zscore = TRUE, min_abs_corr = 0.5)
  }
  compare_guarded_fit <- function(reference_x, reference_y, actual_x, actual_y) {
    set.seed(111L)
    reference <- legacy_fit(reference_x, reference_y)
    reference_rng <- .Random.seed
    set.seed(111L)
    actual <- fitSusieModel(actual_x, actual_y)
    stopifnot("Guarded PIPs must match the legacy fit" =
                max(abs(reference$pip - actual$pip)) < 1e-10,
              "Guarded z-scores must match the legacy fit" =
                max(abs(reference$z - actual$z)) < 1e-10,
              identical(reference$sets$cs, actual$sets$cs),
              identical(reference$converged, actual$converged),
              identical(reference_rng, .Random.seed))
  }
  compare_guarded_fit(legacy_x, guard_hat %*% guard_y, guarded_x,
                      residualizeWithCovariates(guard_adjustment, guard_y))
  near_x <- ordinary_x[, -1L, drop = FALSE]
  compare_guarded_fit(near_hat %*% near_x, legacy_y,
                      residualizeWithCovariates(near_adjustment, near_x), guarded_y)

  set.seed(583L)
  sample_ids <- paste0("s", seq_len(180L))
  covariates <- matrix(rnorm(180L * 3L), nrow = 180L,
                       dimnames = list(sample_ids, paste0("c", seq_len(3L))))
  design <- cbind(covariates, 1)
  values <- matrix(rnorm(180L * 5L), nrow = 180L,
                   dimnames = list(sample_ids, paste0("v", seq_len(5L))))
  dense <- diag(nrow(design)) - design %*% solve(crossprod(design)) %*% t(design)
  adjustment <- createCovariateAdjustment(design)
  residualized <- residualizeWithCovariates(adjustment, values)
  stopifnot(identical(dimnames(residualized), dimnames(dense %*% values)),
            max(abs(residualized - dense %*% values)) < 1e-12)
  vector_values <- values[, 1L]
  vector_result <- residualizeWithCovariates(adjustment, vector_values)
  stopifnot(identical(dim(vector_result), c(180L, 1L)),
            max(abs(vector_result - dense %*% vector_values)) < 1e-12)
  intercept <- matrix(1, nrow = 180L, ncol = 1L)
  stopifnot(max(abs(residualizeWithCovariates(createCovariateAdjustment(intercept), values) -
                     (diag(180L) - matrix(1 / 180L, 180L, 180L)) %*% values)) < 1e-12)
  singular <- cbind(design, design[, 1L])
  stopifnot(inherits(try(createCovariateAdjustment(singular), silent = TRUE), "try-error"))
  large_design <- cbind(rep(1, 10000L), seq_len(10000L), sin(seq_len(10000L)))
  stopifnot(as.numeric(object.size(createCovariateAdjustment(large_design))) < 5 * 1024^2)

  positions <- 1000000L + seq_len(50L) * 10L
  genotype <- matrix(rbinom(50L * 180L, 2L, 0.3), nrow = 50L,
                      dimnames = list(paste0("chr1_", positions, "_A_G"), sample_ids))
  phenotype <- 1.2 * genotype[8L, ] - genotype[32L, ] + rnorm(180L, sd = 0.5)
  phenotype_2 <- genotype[15L, ] + rnorm(180L, sd = 0.5)
  genotype[4L, ] <- 0
  genotype[5L, ] <- -1
  genotype[2L, 3L] <- -1
  genotype[6L, 4L] <- NA_real_
  make_se <- function(positions_for_traits = c(1000200L, 1000250L)) {
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = rbind(gene1 = phenotype, gene2 = phenotype_2)),
      colData = data.frame(sample_id = sample_ids, genotype_id = sample_ids,
                           qtl_group = "ALL", row.names = sample_ids),
      rowData = data.frame(phenotype_id = c("gene1", "gene2"), chromosome = "1",
                           phenotype_pos = positions_for_traits, row.names = c("gene1", "gene2")))
  }
  se <- make_se()
  # Derive expected matrices with the old dense residualization. Do not use
  # the new adjustment or fitting wrapper to build the reference.
  dense_reference <- function(raw_genotypes, trait, samples = sample_ids,
                              covariate_values = covariates,
                              trait_samples = sample_ids,
                              selected_variants = NULL) {
    trait <- qnorm((rank(trait, na.last = "keep") - 0.5) / sum(!is.na(trait)))
    cov <- cbind(covariate_values[samples, , drop = FALSE], 1)
    hat <- diag(length(samples)) - cov %*% solve(crossprod(cov)) %*% t(cov)
    y <- hat %*% trait[match(samples, trait_samples)]
    names(y) <- samples
    gt <- raw_genotypes[, samples, drop = FALSE]
    gt <- gt[rowSums(round(gt, 0), na.rm = TRUE) != 0, , drop = FALSE]
    gt[gt == -1] <- NA_real_
    all_missing <- rowSums(!is.na(gt)) == 0L
    means <- rowMeans(gt, na.rm = TRUE)
    means[is.na(means)] <- 0
    missing <- which(is.na(gt), arr.ind = TRUE)
    gt[missing] <- means[missing[, 1L]]
    frequency <- rowMeans(gt) / 2
    expected_maf <- pmin(frequency, 1 - frequency)
    invalid <- apply(gt, 1L, function(row) any(!is.finite(row) | row < 0 | row > 2))
    expected_maf[all_missing | invalid] <- NA_real_
    x <- hat %*% t(gt - rowMeans(gt))
    if (!is.null(selected_variants)) x <- x[, selected_variants, drop = FALSE]
    fitted <- susieR::susie(x, y, L = 10,
                            estimate_residual_variance = TRUE,
                            estimate_prior_variance = TRUE,
                            scaled_prior_variance = 0.1, verbose = FALSE,
                            compute_univariate_zscore = TRUE, min_abs_corr = 0.5)
    fitted$phenotype_residual_SD <- stats::sd(drop(y))
    fitted$variant_id <- colnames(x)
    fitted$cohort_maf <- unname(expected_maf[match(colnames(x), rownames(gt))])
    fitted
  }
  compare_fit <- function(reference, actual) {
    stopifnot(identical(reference$variant_id, actual$variant_id),
              isTRUE(all.equal(reference$cohort_maf, actual$cohort_maf, tolerance = 1e-12)),
              isTRUE(all.equal(reference$phenotype_residual_SD,
                               actual$phenotype_residual_SD, tolerance = 1e-12)),
              identical(reference$sets$cs, actual$sets$cs),
              identical(reference$converged, actual$converged),
              identical(attributes(reference$z), attributes(actual$z)),
              identical(is.na(reference$z), is.na(actual$z)),
              max(abs(reference$pip - actual$pip)) < 1e-8,
              max(abs(reference$z - actual$z), na.rm = TRUE) < 1e-9)
  }
  covariates_reordered <- covariates[rev(sample_ids), , drop = FALSE]
  set.seed(159L)
  expected <- dense_reference(genotype, phenotype)
  expected_rng <- .Random.seed
  set.seed(159L)
  actual <- finemapPhenotype("gene1", se, genotype, covariates_reordered, 1000L)
  compare_fit(expected, actual)
  stopifnot("An entirely missing legacy predictor must not report observed MAF zero" =
              is.na(actual$cohort_maf[match(rownames(genotype)[[5L]], actual$variant_id)]))
  assert_flag_table(actual$covariate_explained_variants)
  stopifnot(identical(expected_rng, .Random.seed))

  additional_samples <- sample_ids[seq_len(140L)]
  additional <- matrix(rbinom(2L * 140L, 2L, 0.3), nrow = 2L,
                        dimnames = list(c("chr1_1000700_G_T", "chr1_1000800_C_A"), rev(additional_samples)))
  merged <- rbind(genotype[, additional_samples, drop = FALSE], additional[, additional_samples, drop = FALSE])
  set.seed(271L)
  expected_additional <- dense_reference(merged, phenotype, additional_samples)
  set.seed(271L)
  actual_additional <- finemapPhenotype("gene1", se, genotype, covariates_reordered, 1000L,
                                       additional_genotypes = additional)
  compare_fit(expected_additional, actual_additional)

  temp <- tempfile("susier-path-test-")
  dir.create(temp)
  on.exit(unlink(temp, recursive = TRUE), add = TRUE)
  dosage_tsv <- file.path(temp, "dosages.tsv")
  dosage <- data.frame(CHROM = "1", POS = positions, REF = "A", ALT = "G") %>%
    dplyr::bind_cols(as.data.frame(genotype[, rev(sample_ids), drop = FALSE]))
  data.table::fwrite(dosage, dosage_tsv, sep = "\t", quote = FALSE, na = "NA", nThread = 1L)
  dosage_bgz <- Rsamtools::bgzip(dosage_tsv, dest = paste0(dosage_tsv, ".bgz"), overwrite = TRUE)
  Rsamtools::indexTabix(dosage_bgz, seq = 1L, start = 2L, end = 2L, skip = 1L, zeroBased = FALSE)
  reusable <- prepareReusableGenotype(c("gene1", "gene2"), se, dosage_bgz,
                                       covariates_reordered, 1000L)
  stopifnot(!is.null(reusable))
  for (trait in c("gene1", "gene2")) {
    y <- if (trait == "gene1") phenotype else phenotype_2
    set.seed(123L)
    expected_reused <- dense_reference(genotype, y)
    set.seed(123L)
    actual_reused <- finemapPhenotype(trait, se, dosage_bgz, covariates_reordered, 1000L,
                                      reusable_genotype = reusable)
    compare_fit(expected_reused, actual_reused)
  }
  stopifnot(is.null(prepareReusableGenotype(c("gene1", "gene2"), make_se(c(1000200L, 2000200L)),
                                            dosage_bgz, covariates_reordered, 1000L)),
            is.null(prepareReusableGenotype(c("gene1", "gene2"), se, dosage_bgz,
                                            covariates_reordered, 1000L, additional_genotypes = additional)))

  # Exercise coordinated fallback at each production fit boundary. The
  # reusable reference multiplies the complete merged window before subsetting.
  guard_positions <- 1000000L + seq_len(30L) * 50L
  guard_genotype <- t(cbind(guard_raw_x, matrix(rbinom(guard_n * 14L, 2L, 0.3), nrow = guard_n)))
  dimnames(guard_genotype) <- list(paste0("chr1_", guard_positions, "_A_G"), guard_ids)
  guard_covariates <- guard_design[, c("sex", "pc"), drop = FALSE]
  guard_trait_2 <- sin(seq_len(guard_n))
  guard_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = rbind(gene1 = guard_y, gene2 = guard_trait_2)),
    colData = data.frame(sample_id = guard_ids, genotype_id = guard_ids,
                         qtl_group = "ALL", row.names = guard_ids),
    rowData = data.frame(phenotype_id = c("gene1", "gene2"), chromosome = "1",
                         phenotype_pos = c(1000400L, 1001200L), row.names = c("gene1", "gene2")))
  exact_fields <- c("variant_id", "pip", "z", "alpha", "mu", "mu2", "lbf_variable",
                     "sets", "V", "sigma2", "converged", "elbo", "niter",
                     "X_column_scale_factors")
  exact_pipeline_fit <- function(expected, actual) {
    actual_without_metadata <- actual
    actual_without_metadata$covariate_explained_variants <- NULL
    stopifnot("Fallback may add only the requested diagnostics metadata" =
                identical(expected, actual_without_metadata))
    for (field in exact_fields) {
      stopifnot("A legacy fallback fit must retain exact native results" =
                  identical(expected[[field]], actual[[field]]))
    }
  }
  set.seed(191L)
  expected_guard <- dense_reference(guard_genotype, guard_y, guard_ids,
                                    guard_covariates, guard_ids)
  expected_guard_rng <- .Random.seed
  set.seed(191L)
  actual_guard <- finemapPhenotype("gene1", guard_se, guard_genotype,
                                   guard_covariates[rev(guard_ids), , drop = FALSE], 500L)
  exact_pipeline_fit(expected_guard, actual_guard)
  assert_flag_table(actual_guard$covariate_explained_variants, rownames(guard_genotype)[1L])
  stopifnot(identical(expected_guard_rng, .Random.seed))

  guard_additional_ids <- guard_ids[seq_len(140L)]
  guard_additional <- matrix(rbinom(2L * 140L, 2L, 0.3), nrow = 2L,
                              dimnames = list(c("chr1_1001600_C_T", "chr1_1001700_G_A"),
                                              rev(guard_additional_ids)))
  guard_merged <- rbind(guard_genotype[, guard_additional_ids, drop = FALSE],
                        guard_additional[, guard_additional_ids, drop = FALSE])
  set.seed(192L)
  expected_guard_additional <- dense_reference(guard_merged, guard_y, guard_additional_ids,
                                               guard_covariates, guard_ids)
  expected_guard_rng <- .Random.seed
  set.seed(192L)
  actual_guard_additional <- finemapPhenotype("gene1", guard_se, guard_genotype,
                                              guard_covariates, 500L,
                                              additional_genotypes = guard_additional)
  exact_pipeline_fit(expected_guard_additional, actual_guard_additional)
  assert_flag_table(actual_guard_additional$covariate_explained_variants,
                    rownames(guard_genotype)[1L])
  stopifnot(identical(expected_guard_rng, .Random.seed))

  guard_dosage <- data.frame(CHROM = "1", POS = guard_positions, REF = "A", ALT = "G") %>%
    dplyr::bind_cols(as.data.frame(guard_genotype[, rev(guard_ids), drop = FALSE]))
  guard_tsv <- file.path(temp, "guard-dosages.tsv")
  data.table::fwrite(guard_dosage, guard_tsv, sep = "\t", quote = FALSE, nThread = 1L)
  guard_bgz <- Rsamtools::bgzip(guard_tsv, dest = paste0(guard_tsv, ".bgz"), overwrite = TRUE)
  Rsamtools::indexTabix(guard_bgz, seq = 1L, start = 2L, end = 2L, skip = 1L, zeroBased = FALSE)
  selected_guard_variants <- rownames(guard_genotype)[guard_positions <= 1000900L]
  guard_reused <- prepareReusableGenotype(c("gene1", "gene2"), guard_se, guard_bgz,
                                           guard_covariates, 500L)
  set.seed(193L)
  expected_guard_reused <- dense_reference(guard_genotype, guard_y, guard_ids,
                                           guard_covariates, guard_ids, selected_guard_variants)
  expected_guard_rng <- .Random.seed
  set.seed(193L)
  actual_guard_reused <- finemapPhenotype("gene1", guard_se, guard_bgz, guard_covariates,
                                          500L, reusable_genotype = guard_reused)
  exact_pipeline_fit(expected_guard_reused, actual_guard_reused)
  assert_flag_table(actual_guard_reused$covariate_explained_variants,
                    rownames(guard_genotype)[1L])
  stopifnot(identical(expected_guard_rng, .Random.seed))
  guard_reused_other <- finemapPhenotype("gene2", guard_se, guard_bgz, guard_covariates,
                                         500L, reusable_genotype = guard_reused)
  assert_flag_table(guard_reused_other$covariate_explained_variants)

  near_covariates <- near_design[, "pc", drop = FALSE]
  late_reused <- prepareReusableGenotype(c("gene2", "gene1"), guard_se, guard_bgz,
                                          near_covariates, 500L)
  # First fit an ordinary phenotype from this fast cache. A later phenotype
  # must reload the original merged window if its response needs dense math.
  invisible(finemapPhenotype("gene2", guard_se, guard_bgz, near_covariates,
                             500L, reusable_genotype = late_reused))
  set.seed(194L)
  expected_late <- dense_reference(guard_genotype, guard_y, guard_ids,
                                    near_covariates, guard_ids, selected_guard_variants)
  expected_late_rng <- .Random.seed
  set.seed(194L)
  actual_late <- finemapPhenotype("gene1", guard_se, guard_bgz, near_covariates,
                                  500L, reusable_genotype = late_reused)
  exact_pipeline_fit(expected_late, actual_late)
  assert_flag_table(actual_late$covariate_explained_variants)
  stopifnot(identical(expected_late_rng, .Random.seed))
  # The rebuilt result must be cached. Remove the source only inside this
  # temporary fixture; another legacy fit must still use the rebuilt cache.
  stopifnot(file.rename(guard_bgz, paste0(guard_bgz, ".saved")))
  set.seed(194L)
  cached_late <- finemapPhenotype("gene1", guard_se, guard_bgz, near_covariates,
                                  500L, reusable_genotype = late_reused)
  exact_pipeline_fit(expected_late, cached_late)
  assert_flag_table(cached_late$covariate_explained_variants)
  stopifnot(identical(expected_late_rng, .Random.seed))

  # Use a real fit with controlled credible-set membership to check that a
  # nonflagged member inherits its set's flag and overlap handling is unchanged.
  cs_fixture <- actual_guard
  cs_fixture$sets$cs <- list(L1 = c(1L, 2L), L2 = c(3L, 4L), L3 = c(2L, 5L))
  cs_fixture$sets$cs_index <- c(1L, 2L, 3L)
  cs_fixture$sets$purity <- data.frame(min.abs.corr = c(0.8, 0.9, 0.85),
                                      mean.abs.corr = c(0.85, 0.92, 0.9),
                                      median.abs.corr = c(0.85, 0.92, 0.9),
                                      row.names = c("L1", "L2", "L3"))
  annotated <- extractResults(cs_fixture)
  stopifnot(identical(annotated$cs_df$has_covariate_explained_variant, c(TRUE, FALSE)),
            identical(annotated$cs_df$cs_id, c("L1", "L2")),
            all(annotated$variant_df$has_covariate_explained_variant[
              annotated$variant_df$variant_id %in% cs_fixture$variant_id[1:2]]),
            !any(annotated$variant_df$has_covariate_explained_variant[
              !annotated$variant_df$variant_id %in% cs_fixture$variant_id[1:2]]))
  cs_fixture$covariate_explained_variants <- NULL
  old_fit_tables <- extractResults(cs_fixture)
  stopifnot(!any(old_fit_tables$cs_df$has_covariate_explained_variant),
            !any(old_fit_tables$variant_df$has_covariate_explained_variant),
            identical(dplyr::select(annotated$cs_df, -has_covariate_explained_variant),
                      dplyr::select(old_fit_tables$cs_df, -has_covariate_explained_variant)),
            identical(dplyr::select(annotated$variant_df, -has_covariate_explained_variant),
                      dplyr::select(old_fit_tables$variant_df, -has_covariate_explained_variant)),
            identical(annotated$lbf_df, old_fit_tables$lbf_df))

  message("Fine-mapping paths passed: dense reference, normal, additional, reused, and fallback guards")
  invisible(TRUE)
}
