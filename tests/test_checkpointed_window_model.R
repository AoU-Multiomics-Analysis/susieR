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

run_named_test("dosage reader can retain missing cells for mask-specific preparation", {
  dosage <- read_checkpointed_window_dosage(
    fixture_path("window_dosage.tsv"),
    impute = FALSE
  )
  expect_identical_value(
    sum(is.na(dosage$genotype)),
    1L,
    "raw missing dosage count"
  )
  expect_identical_value(
    dosage$genotype[3L, 4L],
    NA_real_,
    "raw missing dosage cell"
  )
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

run_named_test("phenotype reader accepts the real chrom metadata header", {
  phenotypes <- read_checkpointed_window_phenotypes(
    fixture_path("chrom_window_phenotypes.tsv"),
    ordered_manifest
  )
  expect_identical_value(
    phenotypes$metadata$chromosome,
    c("7", "12"),
    "chrom-header chromosome metadata"
  )
  expect_identical_value(
    phenotypes$metadata$start,
    c(700000L, 1200000L),
    "chrom-header start metadata"
  )
  expect_identical_value(
    phenotypes$metadata$end,
    c(700001L, 1200001L),
    "chrom-header end metadata"
  )
  expect_identical_value(
    phenotypes$metadata$phenotype_id,
    c("linked_expression", "linked_splicing"),
    "chrom-header phenotype IDs"
  )
  expect_identical_value(
    unname(phenotypes$values$linked_expression),
    c(-8, -1),
    "chrom-header phenotype values"
  )
})

covariate_paths <- c(
  fixture_path("shared_covariates.tsv"),
  fixture_path("expression_covariates.tsv")
)
covariate_labels <- c("shared", "expression")

run_named_test("covariate reader preserves numeric sample IDs from the header", {
  covariates <- read_checkpointed_covariates(
    fixture_path("numeric_sample_covariates.tsv"),
    "shared"
  )
  expect_identical_value(
    covariates[[1L]]$covariate_ids,
    c("GENETICPC1", "GENETICPC2"),
    "numeric-header covariate IDs"
  )
  expect_identical_value(
    covariates[[1L]]$sample_ids,
    c("1000291", "1000610"),
    "numeric-header sample IDs"
  )
  if (!isTRUE(all.equal(
    unname(covariates[[1L]]$values),
    matrix(c(0.25, 1.5, -0.5, 2.25), nrow = 2L),
    tolerance = 1e-12
  ))) {
    stop("Numeric-header covariate values were not preserved.", call. = FALSE)
  }
})

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

run_named_test("one phenotype NA does not remove a sample from another phenotype", {
  phenotype_data <- data.table::fread(
    fixture_path("window_phenotypes.bed.gz"),
    check.names = FALSE,
    data.table = FALSE
  ) |>
    tibble::as_tibble()
  phenotype_data$sample_05[phenotype_data$phenotype_id == "linked_expression"] <-
    NA_real_
  phenotype_path <- tempfile(fileext = ".tsv")
  on.exit(unlink(phenotype_path), add = TRUE)
  readr::write_tsv(phenotype_data, phenotype_path)
  na_config <- config
  na_config$phenotype_data <- phenotype_path
  na_inputs <- load_checkpointed_window_inputs(na_config, ordered_manifest)

  expect_true_value("sample_05" %in% na_inputs$sample_ids, "shared overlap sample")
  expect_identical_value(
    na_inputs$phenotypes$linked_expression[["sample_05"]],
    NA_real_,
    "phenotype-specific missing value"
  )
  expect_true_value(
    is.finite(na_inputs$phenotypes$linked_splicing[["sample_05"]]),
    "unaffected phenotype sample"
  )
})

run_named_test("sample allowlist requires exactly one named column", {
  invalid_allowlist <- tempfile(fileext = ".tsv")
  on.exit(unlink(invalid_allowlist), add = TRUE)
  readr::write_tsv(
    tibble::tibble(sample_id = c("sample_05", "sample_06"), note = c("a", "b")),
    invalid_allowlist
  )
  expect_error_message(
    read_checkpointed_sample_allowlist(invalid_allowlist),
    "exactly one named"
  )
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

run_named_test("modality covariate gaps affect only that modality cache", {
  splicing_sample_ids <- sprintf("sample_%02d", 6:40)
  splicing_path <- tempfile(fileext = ".tsv")
  on.exit(unlink(splicing_path), add = TRUE)
  splicing_covariates <- tibble::tibble(covariate_id = "splicing_pc1") |>
    dplyr::bind_cols(
      tibble::as_tibble(
        matrix(
          seq_along(splicing_sample_ids) / 10,
          nrow = 1L,
          dimnames = list(NULL, splicing_sample_ids)
        )
      )
    )
  readr::write_tsv(splicing_covariates, splicing_path)
  modality_config <- config
  modality_config$covariate_files <- c(
    modality_config$covariate_files,
    splicing_path
  )
  modality_config$covariate_modalities <- c(
    modality_config$covariate_modalities,
    "splicing"
  )

  modality_inputs <- load_checkpointed_window_inputs(
    modality_config,
    ordered_manifest
  )
  modality_plan <- build_checkpointed_window_cache_plan(
    modality_inputs,
    ordered_manifest
  )
  expression_group <- modality_plan$groups[[
    modality_plan$phenotypes$linked_expression$group_key
  ]]
  splicing_group <- modality_plan$groups[[
    modality_plan$phenotypes$linked_splicing$group_key
  ]]
  expect_true_value(
    "sample_05" %in% modality_inputs$sample_ids,
    "base alignment retains expression sample"
  )
  expect_true_value(
    "sample_05" %in% expression_group$retained_sample_ids,
    "expression cache retains expression sample"
  )
  expect_true_value(
    !"sample_05" %in% splicing_group$retained_sample_ids,
    "splicing cache removes missing-covariate sample"
  )
  expect_identical_value(expression_group$overlap_samples, 36L, "expression overlap")
  expect_identical_value(splicing_group$overlap_samples, 35L, "splicing overlap")

  modality_cache <- new_checkpointed_window_genotype_cache(
    modality_inputs$genotype,
    modality_inputs$variant_info$variant_id,
    modality_plan,
    raw_dosage_samples = modality_inputs$raw_dosage_samples,
    overlap_samples = modality_inputs$overlap_samples
  )
  expression_prepared <- modality_cache$get("linked_expression")
  splicing_prepared <- modality_cache$get("linked_splicing")
  expect_identical_value(expression_prepared$qc$overlap_samples, 36L, "expression QC overlap")
  expect_identical_value(splicing_prepared$qc$overlap_samples, 35L, "splicing QC overlap")
  expect_true_value(
    expression_prepared$cache_key != splicing_prepared$cache_key,
    "modality overlap cache identities"
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

run_named_test("cache key covers samples, retained variants, and covariate design", {
  sample_ids <- model_inputs$sample_ids
  design <- cbind(
    intercept = 1,
    applicable_covariate_matrix(
      model_inputs$covariates,
      "expression",
      sample_ids
    )
  )
  variant_ids <- model_inputs$variant_info$variant_id
  key <- checkpointed_genotype_cache_key(sample_ids, variant_ids, design)
  changed_design <- design
  changed_design[1L, 2L] <- changed_design[1L, 2L] + 1e-6
  changed_keys <- c(
    checkpointed_genotype_cache_key(rev(sample_ids), variant_ids, design[nrow(design):1L, , drop = FALSE]),
    checkpointed_genotype_cache_key(sample_ids, rev(variant_ids), design),
    checkpointed_genotype_cache_key(sample_ids, variant_ids, changed_design)
  )
  expect_true_value(all(changed_keys != key), "cache-key identity inputs")
})

run_named_test("design and cache identities distinguish nearby numeric values", {
  small_value <- 1.1199976554521247e-146
  adjacent_value <- small_value +
    abs(small_value) * 2 * .Machine$double.eps
  first_design <- matrix(
    c(1, small_value, 1, 1),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(NULL, c("intercept", "covariate"))
  )
  adjacent_design <- first_design
  adjacent_design[1L, 2L] <- adjacent_value
  sample_ids <- c("sample_1", "sample_2")
  variant_ids <- "chr1_100_A_G"

  if (identical(
    checkpointed_design_identity(sample_ids, first_design),
    checkpointed_design_identity(sample_ids, adjacent_design)
  )) {
    stop("Distinct designs must have distinct design IDs.", call. = FALSE)
  }
  if (identical(
    checkpointed_genotype_cache_key(sample_ids, variant_ids, first_design),
    checkpointed_genotype_cache_key(sample_ids, variant_ids, adjacent_design)
  )) {
    stop("Distinct designs must have distinct cache keys.", call. = FALSE)
  }
})

run_named_test("same design and mask prepares genotype once", {
  replica_manifest <- dplyr::bind_rows(
    ordered_manifest[1L, ],
    ordered_manifest[1L, ] |>
      dplyr::mutate(
        phenotype_id = "linked_expression_replica",
        processing_index = 1L,
        phenotype_key = checkpoint_phenotype_key(
          "chr1_0_2000000",
          "expression",
          "linked_expression_replica"
        )
      )
  )
  replica_inputs <- model_inputs
  replica_inputs$phenotypes$linked_expression_replica <-
    replica_inputs$phenotypes$linked_expression
  plan <- build_checkpointed_window_cache_plan(replica_inputs, replica_manifest)
  cache <- new_checkpointed_window_genotype_cache(
    replica_inputs$genotype,
    replica_inputs$variant_info$variant_id,
    plan,
    raw_dosage_samples = replica_inputs$raw_dosage_samples,
    overlap_samples = replica_inputs$overlap_samples
  )
  first <- cache$get("linked_expression")
  second <- cache$get("linked_expression_replica")
  metrics <- cache$metrics()

  expect_identical_value(
    plan$phenotypes$linked_expression$group_key,
    plan$phenotypes$linked_expression_replica$group_key,
    "shared preparation key"
  )
  expect_identical_value(metrics$total_preparations, 1L, "genotype preparations")
  expect_identical_value(first$cache_key, second$cache_key, "shared final cache key")
  expect_true_value(metrics$raw_genotype_released, "raw genotype release")
  cache$release("linked_expression")
  expect_identical_value(
    cache$metrics()$resident_entries,
    1L,
    "cache retained before last use"
  )
  cache$release("linked_expression_replica")
  expect_identical_value(
    cache$metrics()$resident_entries,
    0L,
    "cache release after last use"
  )
})

run_named_test("cache planning stores one design for one hundred same-group phenotypes", {
  phenotype_ids <- sprintf("expression_replica_%03d", seq_len(100L))
  repeated_manifest <- purrr::map_dfr(
    seq_along(phenotype_ids),
    function(index) {
      ordered_manifest[1L, ] |>
        dplyr::mutate(
          phenotype_id = phenotype_ids[[index]],
          processing_index = as.integer(index - 1L),
          phenotype_key = checkpoint_phenotype_key(
            "chr1_0_2000000",
            "expression",
            phenotype_ids[[index]]
          )
        )
    }
  )
  repeated_inputs <- model_inputs
  repeated_inputs$phenotypes <- stats::setNames(
    rep(list(model_inputs$phenotypes$linked_expression), 100L),
    phenotype_ids
  )

  original_applicable_covariate_matrix <- applicable_covariate_matrix
  design_preprocessing_count <- 0L
  assign(
    "applicable_covariate_matrix",
    function(...) {
      design_preprocessing_count <<- design_preprocessing_count + 1L
      original_applicable_covariate_matrix(...)
    },
    envir = .GlobalEnv
  )
  on.exit(
    assign(
      "applicable_covariate_matrix",
      original_applicable_covariate_matrix,
      envir = .GlobalEnv
    ),
    add = TRUE
  )

  plan <- build_checkpointed_window_cache_plan(
    repeated_inputs,
    repeated_manifest
  )
  expect_identical_value(
    design_preprocessing_count,
    1L,
    "design preprocessing count"
  )
  expect_identical_value(length(plan$groups), 1L, "stored group count")
  expect_identical_value(length(plan$phenotypes), 100L, "phenotype reference count")
  expect_true_value(
    all(purrr::map_lgl(
      plan$phenotypes,
      ~ identical(names(.x), c(
        "processing_index",
        "phenotype_id",
        "modality",
        "group_key"
      ))
    )),
    "lightweight phenotype references"
  )
  expect_identical_value(
    unique(purrr::map_chr(plan$phenotypes, "group_key")),
    names(plan$groups),
    "single shared group reference"
  )
})

run_named_test("different design or sample mask prepares separately", {
  different_mask_inputs <- model_inputs
  different_mask_inputs$phenotypes$linked_expression[[1L]] <- NA_real_
  mask_plan <- build_checkpointed_window_cache_plan(
    different_mask_inputs,
    ordered_manifest
  )
  mask_cache <- new_checkpointed_window_genotype_cache(
    different_mask_inputs$genotype,
    different_mask_inputs$variant_info$variant_id,
    mask_plan,
    raw_dosage_samples = different_mask_inputs$raw_dosage_samples,
    overlap_samples = different_mask_inputs$overlap_samples
  )
  mask_cache$get("linked_expression")
  expect_true_value(
    !mask_cache$metrics()$raw_genotype_released,
    "raw genotype retained before all preparations"
  )
  mask_cache$get("linked_splicing")
  expect_identical_value(
    mask_cache$metrics()$total_preparations,
    2L,
    "separate genotype preparations"
  )
  expect_true_value(
    mask_plan$phenotypes$linked_expression$group_key !=
      mask_plan$phenotypes$linked_splicing$group_key,
    "different preparation keys"
  )
  expect_true_value(
    mask_cache$metrics()$raw_genotype_released,
    "raw genotype released after all preparations"
  )
})

run_named_test("cached and one-shot scientific results agree", {
  one_record <- ordered_manifest[1L, , drop = FALSE]
  plan <- build_checkpointed_window_cache_plan(model_inputs, one_record)
  cache <- new_checkpointed_window_genotype_cache(
    model_inputs$genotype,
    model_inputs$variant_info$variant_id,
    plan,
    raw_dosage_samples = model_inputs$raw_dosage_samples,
    overlap_samples = model_inputs$overlap_samples
  )
  prepared <- cache$get("linked_expression")
  invisible(utils::capture.output(
    cached <- fit_prepared_checkpointed_window_phenotype(
      prepared,
      model_inputs$phenotypes$linked_expression,
      checkpointed_susie_settings()
    )
  ))
  one_shot_covariates <- applicable_covariate_matrix(
    model_inputs$covariates,
    "expression",
    model_inputs$sample_ids
  )
  invisible(utils::capture.output(
    one_shot <- fit_checkpointed_window_phenotype(
      model_inputs$genotype,
      model_inputs$phenotypes$linked_expression,
      one_shot_covariates,
      model_inputs$variant_info$variant_id,
      checkpointed_susie_settings()
    )
  ))

  expect_identical_value(cached$status, one_shot$status, "fit status")
  expect_identical_value(cached$fit$variant_id, one_shot$fit$variant_id, "fit variants")
  if (!isTRUE(all.equal(cached$fit$pip, one_shot$fit$pip, tolerance = 1e-12))) {
    stop("Cached and one-shot PIPs differ.", call. = FALSE)
  }
  expect_identical_value(cached$qc$raw_dosage_samples, 40L, "raw sample count")
  expect_identical_value(cached$qc$overlap_samples, 36L, "overlap sample count")
  expect_identical_value(cached$qc$phenotype_retained_samples, 36L, "phenotype sample count")
  expect_true_value(nzchar(cached$qc$design_id), "design identity")
})

run_named_test("constant phenotype takes precedence over a genotype skip", {
  prepared <- list(
    skip_reason = "NO_ALTERNATIVE_ALLELE",
    sample_ids = c("sample_a", "sample_b", "sample_c"),
    qc = list(genotype_skip_reason = "NO_ALTERNATIVE_ALLELE")
  )
  constant <- fit_prepared_checkpointed_window_phenotype(
    prepared,
    stats::setNames(rep(2, 3L), prepared$sample_ids),
    checkpointed_susie_settings()
  )
  expect_identical_value(
    constant$reason,
    "ZERO_PHENOTYPE_VARIANCE",
    "constant phenotype exclusion precedence"
  )
  expect_identical_value(
    constant$qc$genotype_skip_reason,
    "NO_ALTERNATIVE_ALLELE",
    "genotype QC preservation"
  )
})

run_named_test("covariate-span phenotype takes precedence over a genotype skip", {
  sample_ids <- sprintf("sample_%02d", seq_len(8L))
  phenotype <- stats::setNames(seq_len(8L), sample_ids)
  rank_normal_scores <- rank_inverse_normal(phenotype)
  design <- cbind(intercept = 1, rank_normal_score = rank_normal_scores)
  rownames(design) <- sample_ids
  prepared <- list(
    skip_reason = "NO_USABLE_VARIANTS",
    sample_ids = sample_ids,
    design_qr = qr(
      design,
      tol = checkpointed_qr_tolerance,
      LAPACK = FALSE
    ),
    qc = list(genotype_skip_reason = "NO_USABLE_VARIANTS")
  )
  result <- fit_prepared_checkpointed_window_phenotype(
    prepared,
    phenotype,
    checkpointed_susie_settings()
  )
  expect_identical_value(
    result$reason,
    "ZERO_PHENOTYPE_VARIANCE",
    "covariate-span phenotype exclusion precedence"
  )
  expect_identical_value(
    result$qc$genotype_skip_reason,
    "NO_USABLE_VARIANTS",
    "covariate-span genotype QC preservation"
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

run_named_test("six-variant credible-set fit uses the fixed ten-effect schema", {
  six_variant_result <- fit_results[[1L]]
  expect_true_value(
    length(six_variant_result$fit$sets$cs) > 0L,
    "six-variant source credible sets"
  )
  tables <- format_checkpointed_susie_tables(
    six_variant_result,
    ordered_manifest[1L, ]
  )
  expect_true_value(nrow(tables$credible_sets) > 0L, "six-variant credible-set rows")
  expect_identical_value(nrow(tables$full_susie), 6L, "six-variant full rows")
  expect_identical_value(nrow(tables$lbf_variable), 6L, "six-variant LBF rows")
  expect_identical_value(
    typeof(tables$full_susie$low_purity),
    "character",
    "nonempty low-purity type"
  )
  expect_true_value(
    all(as.matrix(tables$full_susie[paste0("alpha", 7:10)]) == 0),
    "absent alpha effects"
  )
  expect_true_value(
    all(as.matrix(tables$full_susie[paste0("mu_", 7:10)]) == 0),
    "absent mu effects"
  )
  expect_true_value(
    all(as.matrix(tables$full_susie[paste0("mu2_", 7:10)]) == 0),
    "absent mu2 effects"
  )
  expect_true_value(
    all(is.na(as.matrix(tables$lbf_variable[paste0("lbf_variable", 7:10)]))),
    "absent LBF effects"
  )
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
  expect_identical_value(
    typeof(tables$full_susie$low_purity),
    "character",
    "full-SuSiE low-purity type"
  )
})

run_named_test("result table validator rejects prototype type drift", {
  tables <- checkpointed_empty_susie_tables()
  tables$credible_sets$position <- as.character(tables$credible_sets$position)
  expect_error_message(
    checkpointed_validate_susie_tables(tables),
    "column types"
  )
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
