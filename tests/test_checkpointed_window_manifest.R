source("tests/test_helpers.R")
source("R/utils/CheckpointedWindowSusieFunctions.R")

manifest <- tibble::tribble(
  ~window_id, ~phenotype_id, ~modality, ~phenotype_file, ~p_value,
  "chr1_0_2000000", "trait_b", "splicing", "window.bed.gz", 1e-6,
  "chr1_0_2000000", "trait_a", "expression", "window.bed.gz", 1e-8,
  "chr1_0_2000000", "trait_c", "expression", "window.bed.gz", 1e-6
)

universal_manifest <- tibble::tribble(
  ~window_id, ~outcome_key, ~phenotype_id, ~modality, ~phenotype_file, ~p_value,
  "chr1_0_2000000", "splicing::trait_b", "trait_b", "splicing", "window.bed.gz", 1e-6,
  "chr1_0_2000000", "expression::target", "target", "expression", "window.bed.gz", NA_real_,
  "chr1_0_2000000", "expression::trait_a", "trait_a", "expression", "window.bed.gz", 1e-8
)

run_named_test("manifest is ordered by p value, modality, and phenotype", {
  ordered <- validate_window_phenotype_manifest(manifest, "chr1_0_2000000")
  expect_identical_value(ordered$phenotype_id, c("trait_a", "trait_c", "trait_b"), "phenotype order")
  expect_identical_value(ordered$processing_index, 0:2, "processing indices")
  expect_identical_value(length(unique(ordered$phenotype_key)), 3L, "phenotype key count")
})

run_named_test("universal manifest selects linked rows and uses outcome keys", {
  ordered <- validate_window_phenotype_manifest(
    universal_manifest,
    "chr1_0_2000000"
  )
  expect_identical_value(
    ordered$phenotype_id,
    c("expression::trait_a", "splicing::trait_b"),
    "universal phenotype order"
  )
  expect_identical_value(ordered$processing_index, 0:1, "processing indices")
})

run_named_test("universal manifest rejects outcome keys that do not match", {
  invalid_key <- dplyr::mutate(
    universal_manifest,
    outcome_key = dplyr::if_else(
      .data$phenotype_id == "trait_a",
      "expression::wrong",
      .data$outcome_key
    )
  )
  expect_error_message(
    validate_window_phenotype_manifest(invalid_key, "chr1_0_2000000"),
    "outcome_key must equal modality::phenotype_id"
  )
})

run_named_test("universal manifest requires at least one linked phenotype", {
  target_only <- dplyr::mutate(universal_manifest, p_value = NA_real_)
  expect_error_message(
    validate_window_phenotype_manifest(target_only, "chr1_0_2000000"),
    "does not contain linked phenotypes with finite p_value values"
  )
})

run_named_test("universal manifest rejects nonfinite and nonnumeric P values", {
  invalid_values <- c("NaN", "Inf", "-Inf", "not-a-number", "-0.1", "1.1")
  purrr::walk(invalid_values, function(invalid_value) {
    invalid_manifest <- dplyr::mutate(
      universal_manifest,
      p_value = as.character(.data$p_value)
    )
    invalid_manifest$p_value[[1L]] <- invalid_value
    expect_error_message(
      validate_window_phenotype_manifest(
        invalid_manifest,
        "chr1_0_2000000"
      ),
      "p_value values must be finite and between zero and one"
    )
  })
})

run_named_test("universal manifest accepts P-value boundaries", {
  boundary_manifest <- universal_manifest
  boundary_manifest$p_value[c(1L, 3L)] <- c(1, 0)
  ordered <- validate_window_phenotype_manifest(
    boundary_manifest,
    "chr1_0_2000000"
  )
  expect_identical_value(ordered$p_value, c(0, 1), "P-value boundaries")
})

run_named_test("outcome keys disambiguate shared IDs across modalities", {
  collision_manifest <- tibble::tribble(
    ~window_id, ~outcome_key, ~phenotype_id, ~modality, ~phenotype_file, ~p_value,
    "chr1_0_2000000", "expression::shared", "shared", "expression", "window.tsv", 1e-8,
    "chr1_0_2000000", "splicing::shared", "shared", "splicing", "window.tsv", 1e-7
  )
  phenotype_path <- tempfile(fileext = ".tsv")
  on.exit(unlink(phenotype_path), add = TRUE)
  readr::write_tsv(
    tibble::tribble(
      ~chrom, ~start, ~end, ~phenotype_id, ~sample_1,
      "chr1", 1L, 2L, "expression::shared", 0.1,
      "chr1", 3L, 4L, "splicing::shared", 0.2
    ),
    phenotype_path
  )
  ordered <- validate_window_phenotype_manifest(
    collision_manifest,
    "chr1_0_2000000"
  )
  phenotypes <- read_checkpointed_window_phenotypes(
    phenotype_path,
    ordered
  )
  expect_identical_value(
    names(phenotypes$values),
    c("expression::shared", "splicing::shared"),
    "cross-modality phenotype IDs"
  )
  expect_identical_value(
    ordered$modality,
    c("expression", "splicing"),
    "cross-modality covariate labels"
  )
})

run_named_test("manifest rejects missing required columns", {
  expect_error_message(
    validate_window_phenotype_manifest(dplyr::select(manifest, -p_value), "chr1_0_2000000"),
    "missing required columns: p_value"
  )
})

run_named_test("manifest rejects duplicate phenotype identifiers", {
  expect_error_message(
    validate_window_phenotype_manifest(dplyr::bind_rows(manifest, manifest[1, ]), "chr1_0_2000000"),
    "duplicate phenotype_id"
  )
})

run_named_test("manifest rejects invalid P values", {
  invalid_p_values <- dplyr::bind_rows(
    dplyr::mutate(manifest[1, ], p_value = -0.1),
    dplyr::mutate(manifest[2, ], p_value = 1.1),
    dplyr::mutate(manifest[3, ], p_value = NA_real_)
  )
  expect_error_message(
    validate_window_phenotype_manifest(invalid_p_values, "chr1_0_2000000"),
    "p_value values must be finite and between zero and one"
  )
})

run_named_test("manifest rejects multiple window identifiers", {
  multiple_windows <- dplyr::mutate(manifest, window_id = c("chr1_0_2000000", "chr2_0_2000000", "chr1_0_2000000"))
  expect_error_message(
    validate_window_phenotype_manifest(multiple_windows, "chr1_0_2000000"),
    "window_id does not match the requested window"
  )
})

run_named_test("manifest rejects empty identifiers", {
  empty_identifier <- dplyr::mutate(manifest, phenotype_id = c("trait_a", "", "trait_c"))
  expect_error_message(
    validate_window_phenotype_manifest(empty_identifier, "chr1_0_2000000"),
    "identifiers cannot be empty"
  )
})

run_named_test("manifest rejects requested-window mismatch", {
  expect_error_message(
    validate_window_phenotype_manifest(manifest, "chr1_0_3000000"),
    "window_id does not match the requested window"
  )
})

run_named_test("settings and phenotype keys are deterministic", {
  settings <- checkpointed_susie_settings()
  expect_identical_value(settings$L, 10L, "L setting")
  expect_identical_value(settings$estimate_residual_variance, TRUE, "residual variance setting")
  expect_identical_value(settings$estimate_prior_variance, TRUE, "prior variance setting")
  expect_identical_value(settings$scaled_prior_variance, 0.1, "scaled prior variance setting")
  expect_identical_value(settings$compute_univariate_zscore, TRUE, "univariate z-score setting")
  expect_identical_value(settings$min_abs_corr, 0.5, "minimum correlation setting")
  expect_identical_value(
    checkpoint_phenotype_key("chr1_0_2000000", "expression", "trait_a"),
    "bf3478b7965f3143af373e20df5e92f336254263cf398801ed34b9cf48633982",
    "phenotype key"
  )
})

run_named_test("file hash is deterministic", {
  test_file <- tempfile(fileext = ".txt")
  on.exit(unlink(test_file), add = TRUE)
  writeBin(charToRaw("canonical content"), test_file)
  expect_identical_value(
    sha256_file(test_file),
    "1ef9517fb8d602685149aea524f214becb9107adc1f79a17e4605f000e823b2c",
    "file hash"
  )
})

run_named_test("analysis ID is canonical and sensitive to identity inputs", {
  ordered <- validate_window_phenotype_manifest(manifest, "chr1_0_2000000")
  settings <- checkpointed_susie_settings()
  input_hashes <- c(genotype = "genotype-hash", phenotype = "phenotype-hash")
  id <- build_checkpoint_analysis_id(
    input_hashes,
    ordered,
    settings,
    container_digest = "sha256:container",
    wrapper_hash = "wrapper-hash"
  )
  expect_identical_value(
    id,
    build_checkpoint_analysis_id(
      input_hashes[c("phenotype", "genotype")],
      ordered,
      settings[rev(names(settings))],
      container_digest = "sha256:container",
      wrapper_hash = "wrapper-hash"
    ),
    "canonical analysis ID"
  )

  changed_p_value <- dplyr::mutate(ordered, p_value = ifelse(.data$phenotype_id == "trait_a", 2e-8, .data$p_value))
  changed_file_hash <- input_hashes
  changed_file_hash[["phenotype"]] <- "changed-phenotype-hash"
  changed_settings <- settings
  changed_settings$L <- 11L
  changed_ids <- c(
    build_checkpoint_analysis_id(input_hashes, changed_p_value, settings, "sha256:container", "wrapper-hash"),
    build_checkpoint_analysis_id(changed_file_hash, ordered, settings, "sha256:container", "wrapper-hash"),
    build_checkpoint_analysis_id(input_hashes, ordered, changed_settings, "sha256:container", "wrapper-hash"),
    build_checkpoint_analysis_id(input_hashes, ordered, settings, "sha256:changed-container", "wrapper-hash")
  )
  if (any(changed_ids == id)) {
    stop("Each changed analysis identity input must change the analysis ID.", call. = FALSE)
  }
})

run_named_test("analysis ID distinguishes nearby very small P values", {
  ordered <- validate_window_phenotype_manifest(manifest, "chr1_0_2000000")
  small_p_value <- 1.1199976554521247e-146
  adjacent_p_value <- small_p_value +
    abs(small_p_value) * 2 * .Machine$double.eps
  first_manifest <- ordered |>
    dplyr::mutate(
      p_value = dplyr::if_else(
        .data$phenotype_id == "trait_a",
        small_p_value,
        .data$p_value
      )
    )
  adjacent_manifest <- first_manifest |>
    dplyr::mutate(
      p_value = dplyr::if_else(
        .data$phenotype_id == "trait_a",
        adjacent_p_value,
        .data$p_value
      )
    )
  identity_arguments <- list(
    input_hashes = c(genotype = "genotype-hash"),
    settings = checkpointed_susie_settings(),
    container_digest = "sha256:container",
    wrapper_hash = "wrapper-hash"
  )
  first_id <- build_checkpoint_analysis_id(
    identity_arguments$input_hashes,
    first_manifest,
    identity_arguments$settings,
    identity_arguments$container_digest,
    identity_arguments$wrapper_hash
  )
  adjacent_id <- build_checkpoint_analysis_id(
    identity_arguments$input_hashes,
    adjacent_manifest,
    identity_arguments$settings,
    identity_arguments$container_digest,
    identity_arguments$wrapper_hash
  )
  if (identical(first_id, adjacent_id)) {
    stop("Distinct P values must produce distinct analysis IDs.", call. = FALSE)
  }
})

run_named_test("analysis ID covers every mutable runtime source", {
  ordered <- validate_window_phenotype_manifest(manifest, "chr1_0_2000000")
  source_hashes <- c(
    runner_wrapper = paste(rep("a", 64L), collapse = ""),
    checkpointed_functions = paste(rep("b", 64L), collapse = ""),
    checkpoint_store = paste(rep("c", 64L), collapse = "")
  )
  analysis_id <- build_checkpoint_analysis_id(
    input_hashes = c(genotype = "genotype-hash"),
    ordered_manifest = ordered,
    settings = checkpointed_susie_settings(),
    container_digest = "sha256:base",
    wrapper_hash = source_hashes[["runner_wrapper"]],
    source_hashes = source_hashes,
    runtime_image = "runner@example@sha256:published"
  )

  changed_functions <- source_hashes
  changed_functions[["checkpointed_functions"]] <- paste(rep("d", 64L), collapse = "")
  changed_store <- source_hashes
  changed_store[["checkpoint_store"]] <- paste(rep("e", 64L), collapse = "")
  changed_ids <- c(
    build_checkpoint_analysis_id(
      c(genotype = "genotype-hash"),
      ordered,
      checkpointed_susie_settings(),
      "sha256:base",
      source_hashes[["runner_wrapper"]],
      source_hashes = changed_functions,
      runtime_image = "runner@example@sha256:published"
    ),
    build_checkpoint_analysis_id(
      c(genotype = "genotype-hash"),
      ordered,
      checkpointed_susie_settings(),
      "sha256:base",
      source_hashes[["runner_wrapper"]],
      source_hashes = changed_store,
      runtime_image = "runner@example@sha256:published"
    )
  )
  expect_true_value <- function(value, label) {
    if (!isTRUE(value)) {
      stop("Expected true for ", label, ".", call. = FALSE)
    }
  }
  expect_true_value(
    all(changed_ids != analysis_id),
    "helper hashes changing the analysis ID"
  )
})

run_named_test("manifest phenotype file matches the configured phenotype data", {
  expect_error_message(
    validate_window_phenotype_manifest(
      dplyr::mutate(manifest, phenotype_file = "other.bed.gz"),
      "chr1_0_2000000",
      phenotype_data = "/inputs/window.bed.gz"
    ),
    "phenotype_file must match"
  )
  expect_error_message(
    validate_window_phenotype_manifest(
      dplyr::mutate(
        manifest,
        phenotype_file = c("window.bed.gz", "other.bed.gz", "window.bed.gz")
      ),
      "chr1_0_2000000",
      phenotype_data = "/inputs/window.bed.gz"
    ),
    "one phenotype_file"
  )
})
