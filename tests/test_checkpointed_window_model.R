source("tests/test_helpers.R")
source("R/utils/InitFunctions.R")
source("R/utils/SusieFunctions.R")
source("R/utils/CheckpointedWindowSusieFunctions.R")

expect_true_value <- function(value, label) {
  if (!isTRUE(value)) {
    stop("Expected true for ", label, ".", call. = FALSE)
  }
  invisible(TRUE)
}

fixture_path <- function(name) {
  file.path("tests", "fixtures", "checkpointed_window", name)
}

raw_manifest <- readr::read_tsv(
  fixture_path("window_phenotypes.tsv"),
  show_col_types = FALSE
)
ordered_manifest <- validate_window_phenotype_manifest(
  raw_manifest,
  "chr1_0_2000000"
)

run_named_test("dosage reader imputes missing values and retains window order", {
  dosage <- read_checkpointed_window_dosage(fixture_path("window_dosage.tsv"))
  expect_identical_value(nrow(dosage$genotype), 6L, "variant count")
  expect_true_value(all(is.finite(dosage$genotype)), "finite imputed dosage")
  expect_identical_value(
    dosage$variant_info$variant_id,
    c(
      "1_100000_A_G", "1_200000_C_T", "1_300000_G_A",
      "1_400000_T_C", "1_500000_A_C", "1_600000_C_G"
    ),
    "prepared-window variant order"
  )
  if (!isTRUE(all.equal(dosage$genotype[3L, 4L], 31 / 39, tolerance = 1e-12))) {
    stop("The -1 dosage cell did not receive the hand-calculated row mean.", call. = FALSE)
  }
})

run_named_test("dosage reader rejects finite values outside zero through two", {
  invalid_path <- tempfile(fileext = ".tsv")
  on.exit(unlink(invalid_path), add = TRUE)
  invalid_dosage <- readr::read_tsv(
    fixture_path("window_dosage.tsv"),
    show_col_types = FALSE
  )
  invalid_dosage$sample_01[[1L]] <- 2.1
  readr::write_tsv(invalid_dosage, invalid_path)
  expect_error_message(
    read_checkpointed_window_dosage(invalid_path),
    "between zero and two"
  )
})

run_named_test("phenotype reader follows the authoritative manifest", {
  phenotypes <- read_checkpointed_window_phenotypes(
    fixture_path("window_phenotypes.bed.gz"),
    ordered_manifest
  )
  expect_identical_value(
    names(phenotypes$values),
    c("linked_expression", "linked_splicing"),
    "phenotype value order"
  )
  expect_identical_value(
    phenotypes$metadata$chromosome,
    c("7", "12"),
    "trans phenotype chromosome order"
  )
  expect_identical_value(length(phenotypes$values$linked_expression), 40L, "phenotype sample count")
})

covariate_paths <- c(
  fixture_path("shared_covariates.tsv"),
  fixture_path("expression_covariates.tsv")
)
covariate_labels <- c("shared", "expression")

run_named_test("covariate selection is modality specific and ordered", {
  covariates <- read_checkpointed_covariates(covariate_paths, covariate_labels)
  sample_ids <- sprintf("sample_%02d", seq_len(40L))
  expression_matrix <- applicable_covariate_matrix(
    covariates,
    "expression",
    sample_ids
  )
  splicing_matrix <- applicable_covariate_matrix(
    covariates,
    "splicing",
    sample_ids
  )
  expect_identical_value(
    colnames(expression_matrix),
    c("shared_pc1", "expression_pc1"),
    "expression covariate order"
  )
  expect_identical_value(colnames(splicing_matrix), "shared_pc1", "splicing covariates")
  if (!isTRUE(all.equal(
    unname(expression_matrix[1L, ]),
    c(-1, sin(1 / 4)),
    tolerance = 1e-12
  ))) {
    stop("Covariate values do not match the deterministic fixture.", call. = FALSE)
  }
})

run_named_test("covariate selection removes constant and dependent columns", {
  sample_ids <- sprintf("sample_%02d", seq_len(40L))
  shared <- readr::read_tsv(covariate_paths[[1L]], show_col_types = FALSE)
  rank_test <- tibble::tibble(
    covariate_id = c("constant_pc", "dependent_pc")
  ) |>
    dplyr::bind_cols(
      tibble::as_tibble(
        rbind(rep(3, length(sample_ids)), 2 * as.numeric(shared[1L, -1L])),
        .name_repair = ~ sample_ids
      )
    )
  rank_path <- tempfile(fileext = ".tsv")
  on.exit(unlink(rank_path), add = TRUE)
  readr::write_tsv(rank_test, rank_path)
  covariates <- read_checkpointed_covariates(
    c(covariate_paths[[1L]], rank_path),
    c("shared", "expression")
  )
  matrix <- applicable_covariate_matrix(covariates, "expression", sample_ids)
  expect_identical_value(colnames(matrix), "shared_pc1", "full-rank basis")
  expect_identical_value(
    attr(matrix, "removed_covariates"),
    c("constant_pc", "dependent_pc"),
    "removed covariates"
  )
})

run_named_test("covariate selection includes the intercept in rank reduction", {
  sample_ids <- sprintf("sample_%02d", seq_len(40L))
  x <- seq(0, 1, length.out = length(sample_ids))
  complement_path <- tempfile(fileext = ".tsv")
  on.exit(unlink(complement_path), add = TRUE)
  complement_covariates <- tibble::tibble(
    covariate_id = c("x", "one_minus_x")
  ) |>
    dplyr::bind_cols(
      tibble::as_tibble(
        rbind(x, 1 - x),
        .name_repair = ~ sample_ids
      )
    )
  readr::write_tsv(complement_covariates, complement_path)
  covariates <- read_checkpointed_covariates(complement_path, "shared")
  matrix <- applicable_covariate_matrix(covariates, "expression", sample_ids)
  expect_identical_value(colnames(matrix), "x", "intercept-aware basis")
  expect_identical_value(
    attr(matrix, "removed_covariates"),
    "one_minus_x",
    "intercept-dependent covariate"
  )
})

run_named_test("covariate selection rejects duplicate covariate IDs", {
  duplicate_path <- tempfile(fileext = ".tsv")
  on.exit(unlink(duplicate_path), add = TRUE)
  duplicate <- readr::read_tsv(covariate_paths[[1L]], show_col_types = FALSE)
  readr::write_tsv(duplicate, duplicate_path)
  covariates <- read_checkpointed_covariates(
    c(covariate_paths[[1L]], duplicate_path),
    c("shared", "expression")
  )
  expect_error_message(
    applicable_covariate_matrix(
      covariates,
      "expression",
      sprintf("sample_%02d", seq_len(40L))
    ),
    "duplicate covariate IDs"
  )
})

config <- list(
  window_dosage = fixture_path("window_dosage.tsv"),
  phenotype_data = fixture_path("window_phenotypes.bed.gz"),
  covariate_files = covariate_paths,
  covariate_modalities = covariate_labels,
  keep_samples = fixture_path("keep_samples.tsv")
)
model_inputs <- load_checkpointed_window_inputs(config, ordered_manifest)

run_named_test("input loader aligns every input in prepared dosage order", {
  expected_samples <- sprintf("sample_%02d", 5:40)
  expect_identical_value(model_inputs$sample_ids, expected_samples, "aligned sample order")
  expect_identical_value(colnames(model_inputs$genotype), expected_samples, "genotype samples")
  expect_identical_value(
    names(model_inputs$phenotypes$linked_expression),
    expected_samples,
    "phenotype samples"
  )
  expect_identical_value(dim(model_inputs$genotype), c(6L, 36L), "aligned genotype dimensions")
})

run_named_test("input loader rejects a manifest modality without covariates", {
  incomplete_config <- config
  incomplete_config$covariate_files <- covariate_paths[[2L]]
  incomplete_config$covariate_modalities <- "expression"
  expect_error_message(
    load_checkpointed_window_inputs(incomplete_config, ordered_manifest),
    "no applicable covariate input"
  )
})

run_named_test("input loader ignores samples from inapplicable covariate files", {
  unrelated_path <- tempfile(fileext = ".tsv")
  on.exit(unlink(unrelated_path), add = TRUE)
  readr::write_tsv(
    tibble::tibble(covariate_id = "unused_pc", sample_05 = 1),
    unrelated_path
  )
  extra_config <- config
  extra_config$covariate_files <- c(extra_config$covariate_files, unrelated_path)
  extra_config$covariate_modalities <- c(
    extra_config$covariate_modalities,
    "proteomics"
  )
  extra_inputs <- load_checkpointed_window_inputs(extra_config, ordered_manifest)
  expect_identical_value(
    extra_inputs$sample_ids,
    sprintf("sample_%02d", 5:40),
    "samples after inapplicable covariate"
  )
})

fit_results <- purrr::map(seq_len(nrow(ordered_manifest)), function(index) {
  record <- ordered_manifest[index, ]
  fit_checkpointed_window_phenotype(
    genotype = model_inputs$genotype,
    phenotype = model_inputs$phenotypes[[record$phenotype_id[[1L]]]],
    covariates = applicable_covariate_matrix(
      model_inputs$covariates,
      record$modality[[1L]],
      model_inputs$sample_ids
    ),
    variant_ids = model_inputs$variant_info$variant_id,
    settings = checkpointed_susie_settings()
  )
})

run_named_test("trans-linked phenotypes use the full usable window", {
  dosage_chromosomes <- unique(model_inputs$variant_info$CHROM)
  phenotype_chromosomes <- unique(model_inputs$phenotype_metadata$chromosome)
  expect_true_value(
    all(!phenotype_chromosomes %in% dosage_chromosomes),
    "trans phenotype coordinates"
  )
  purrr::walk(fit_results, function(result) {
    expect_true_value(result$status %in% c("CONVERGED", "NONCONVERGED"), "fit status")
    expect_identical_value(
      result$fit$variant_id,
      model_inputs$variant_info$variant_id,
      "full-window variant order"
    )
    expect_true_value(validate_checkpointed_susie_fit(result$fit), "fit structure")
    expect_true_value(is.vector(result$fit$pip), "vector PIP")
    expect_true_value(nrow(result$fit$alpha) > 0L, "nonempty alpha")
    expect_true_value(nrow(result$fit$mu) > 0L, "nonempty mu")
    expect_true_value(nrow(result$fit$mu2) > 0L, "nonempty mu2")
    expect_true_value(
      is.logical(result$fit$converged) && length(result$fit$converged) == 1L,
      "scalar logical convergence"
    )
    expect_identical_value(result$qc$input_variants, 6L, "input variants")
    expect_identical_value(result$qc$retained_variants, 6L, "retained variants")
  })
  expression_pip <- fit_results[[1L]]$fit$pip
  expect_identical_value(
    sort(order(expression_pip, decreasing = TRUE)[1:2]),
    c(2L, 5L),
    "top expression-trait variants"
  )
})

empty_covariates <- matrix(
  numeric(36L * 0L),
  nrow = 36L,
  dimnames = list(model_inputs$sample_ids, character())
)

run_named_test("fit adapter rejects a phenotype in the covariate span", {
  sample_count <- length(model_inputs$sample_ids)
  phenotype <- seq_len(sample_count)
  rank_normal_scores <- stats::qnorm(
    (seq_len(sample_count) - 0.5) / sample_count
  )
  covariates <- matrix(
    rank_normal_scores,
    ncol = 1L,
    dimnames = list(model_inputs$sample_ids, "rank_normal_score")
  )
  result <- fit_checkpointed_window_phenotype(
    model_inputs$genotype,
    phenotype,
    covariates,
    model_inputs$variant_info$variant_id,
    checkpointed_susie_settings()
  )
  expect_identical_value(result$status, "SKIPPED", "span-phenotype status")
  expect_identical_value(
    result$reason,
    "ZERO_PHENOTYPE_VARIANCE",
    "span-phenotype reason"
  )
})

run_named_test("fit adapter rejects variants in the covariate span", {
  sample_count <- length(model_inputs$sample_ids)
  x <- seq(0, 1, length.out = sample_count)
  covariates <- matrix(
    x,
    ncol = 1L,
    dimnames = list(model_inputs$sample_ids, "x")
  )
  genotype <- rbind(
    in_span_1 = 0.25 + x,
    in_span_2 = 1.75 - x
  )
  result <- fit_checkpointed_window_phenotype(
    genotype,
    model_inputs$phenotypes$linked_expression,
    covariates,
    rownames(genotype),
    checkpointed_susie_settings()
  )
  expect_identical_value(result$status, "SKIPPED", "span-variant status")
  expect_identical_value(
    result$reason,
    "NO_USABLE_VARIANTS",
    "span-variant reason"
  )
  expect_identical_value(result$qc$retained_variants, 0L, "span variants retained")
  expect_identical_value(
    result$qc$removed_variant_ids,
    c("in_span_1", "in_span_2"),
    "span variants removed"
  )
})

run_named_test("fit adapter returns exact expected exclusion reasons", {
  zero_variance <- fit_checkpointed_window_phenotype(
    model_inputs$genotype,
    rep(2, 36L),
    empty_covariates,
    model_inputs$variant_info$variant_id,
    checkpointed_susie_settings()
  )
  expect_identical_value(zero_variance$status, "SKIPPED", "zero-variance status")
  expect_identical_value(zero_variance$reason, "ZERO_PHENOTYPE_VARIANCE", "zero-variance reason")

  no_variable <- fit_checkpointed_window_phenotype(
    matrix(1, nrow = 6L, ncol = 36L),
    model_inputs$phenotypes$linked_expression,
    empty_covariates,
    model_inputs$variant_info$variant_id,
    checkpointed_susie_settings()
  )
  expect_identical_value(no_variable$reason, "NO_USABLE_VARIANTS", "no-variable reason")
  expect_identical_value(no_variable$status, "SKIPPED", "no-variable status")

  no_alternative <- fit_checkpointed_window_phenotype(
    matrix(0, nrow = 6L, ncol = 36L),
    model_inputs$phenotypes$linked_expression,
    empty_covariates,
    model_inputs$variant_info$variant_id,
    checkpointed_susie_settings()
  )
  expect_identical_value(no_alternative$reason, "NO_ALTERNATIVE_ALLELE", "no-alternative reason")
  expect_identical_value(no_alternative$status, "SKIPPED", "no-alternative status")

  no_usable <- fit_checkpointed_window_phenotype(
    matrix(NA_real_, nrow = 6L, ncol = 36L),
    model_inputs$phenotypes$linked_expression,
    empty_covariates,
    model_inputs$variant_info$variant_id,
    checkpointed_susie_settings()
  )
  expect_identical_value(no_usable$reason, "NO_USABLE_VARIANTS", "no-usable reason")
  expect_identical_value(no_usable$status, "SKIPPED", "no-usable status")

  too_few <- fit_checkpointed_window_phenotype(
    model_inputs$genotype[, 1L, drop = FALSE],
    model_inputs$phenotypes$linked_expression[[1L]],
    matrix(numeric(), nrow = 1L, ncol = 0L),
    model_inputs$variant_info$variant_id,
    checkpointed_susie_settings()
  )
  expect_identical_value(too_few$reason, "TOO_FEW_ALIGNED_SAMPLES", "sample-count reason")
  expect_identical_value(too_few$status, "SKIPPED", "sample-count status")
})

credible_set_columns <- c(
  "molecular_trait_id", "variant", "chromosome", "position", "ref", "alt",
  "cs_id", "cs_index", "region", "pip", "z", "cs_min_r2", "cs_avg_r2",
  "cs_size", "posterior_mean", "posterior_sd", "cs_log10bf"
)
lbf_columns <- c(
  "molecular_trait_id", "region", "variant", "chromosome", "position",
  paste0("lbf_variable", seq_len(10L))
)
full_susie_columns <- c(
  "molecular_trait_id", "variant", "chromosome", "position", "ref", "alt",
  "cs_id", "cs_index", "low_purity", "region", "pip", "z",
  "posterior_mean", "posterior_sd", "X_column_scale_factors",
  paste0("alpha", seq_len(10L)), paste0("mu_", seq_len(10L)),
  paste0("mu2_", seq_len(10L))
)

run_named_test("formatter preserves all three result schemas", {
  z_dimensions_before <- dim(fit_results[[1L]]$fit$z)
  tables <- format_checkpointed_susie_tables(fit_results[[1L]], ordered_manifest[1L, ])
  expect_identical_value(
    dim(fit_results[[1L]]$fit$z),
    z_dimensions_before,
    "saved fit z dimensions"
  )
  expect_identical_value(names(tables), c("credible_sets", "lbf_variable", "full_susie"), "table names")
  expect_identical_value(names(tables$credible_sets), credible_set_columns, "credible-set columns")
  expect_identical_value(names(tables$lbf_variable), lbf_columns, "LBF columns")
  expect_identical_value(names(tables$full_susie), full_susie_columns, "full-SuSiE columns")
  if (nrow(tables$full_susie) > 0L) {
    expect_true_value(
      all(tables$full_susie$molecular_trait_id == "linked_expression"),
      "formatted phenotype ID"
    )
    expect_true_value(
      all(tables$full_susie$region == "chr1_0_2000000"),
      "formatted window ID"
    )
    expect_identical_value(
      tables$full_susie$variant,
      model_inputs$variant_info$variant_id,
      "formatted full-window variant order"
    )
  }

  output_dir <- tempfile("checkpointed-model-tables-")
  paths <- write_checkpointed_susie_tables(tables, output_dir)
  expect_true_value(all(file.exists(unlist(paths))), "written result tables")
  written <- purrr::map(paths, arrow::read_parquet)
  expect_identical_value(names(written$credible_sets), credible_set_columns, "written credible-set columns")
  expect_identical_value(names(written$lbf_variable), lbf_columns, "written LBF columns")
  expect_identical_value(names(written$full_susie), full_susie_columns, "written full-SuSiE columns")
})

run_named_test("skipped formatter writes typed empty schemas", {
  skipped <- list(
    status = "SKIPPED",
    fit = NULL,
    reason = "NO_USABLE_VARIANTS",
    qc = list()
  )
  tables <- format_checkpointed_susie_tables(skipped, ordered_manifest[1L, ])
  expect_identical_value(
    unname(vapply(tables, nrow, integer(1L))),
    c(0L, 0L, 0L),
    "empty row counts"
  )
  expect_identical_value(names(tables$credible_sets), credible_set_columns, "empty credible-set columns")
  expect_identical_value(names(tables$lbf_variable), lbf_columns, "empty LBF columns")
  expect_identical_value(names(tables$full_susie), full_susie_columns, "empty full-SuSiE columns")
  expect_identical_value(typeof(tables$credible_sets$position), "integer", "credible-set position type")
  expect_identical_value(typeof(tables$lbf_variable$position), "integer", "LBF position type")
  expect_identical_value(typeof(tables$full_susie$pip), "double", "full-SuSiE PIP type")
})

run_named_test("valid fit without credible sets writes typed empty schemas", {
  fit_without_sets <- fit_results[[1L]]$fit
  fit_without_sets$sets$cs <- list()
  fit_without_sets$sets$cs_index <- integer()
  fit_without_sets$sets$purity <- matrix(
    numeric(),
    nrow = 0L,
    ncol = 3L,
    dimnames = list(
      character(),
      c("min.abs.corr", "mean.abs.corr", "median.abs.corr")
    )
  )
  result_without_sets <- fit_results[[1L]]
  result_without_sets$fit <- fit_without_sets
  tables <- format_checkpointed_susie_tables(
    result_without_sets,
    ordered_manifest[1L, ]
  )
  expect_identical_value(
    unname(vapply(tables, nrow, integer(1L))),
    c(0L, 0L, 0L),
    "no-credible-set row counts"
  )
  expect_identical_value(names(tables$credible_sets), credible_set_columns, "no-set credible schema")
  expect_identical_value(names(tables$lbf_variable), lbf_columns, "no-set LBF schema")
  expect_identical_value(names(tables$full_susie), full_susie_columns, "no-set full schema")
})
