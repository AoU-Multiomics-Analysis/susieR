source("tests/test_helpers.R")
source("R/utils/CheckpointedWindowSusieFunctions.R")

manifest <- tibble::tribble(
  ~window_id, ~phenotype_id, ~modality, ~phenotype_file, ~p_value,
  "chr1_0_2000000", "trait_b", "splicing", "window.bed.gz", 1e-6,
  "chr1_0_2000000", "trait_a", "expression", "window.bed.gz", 1e-8,
  "chr1_0_2000000", "trait_c", "expression", "window.bed.gz", 1e-6
)

run_named_test("manifest is ordered by p value, modality, and phenotype", {
  ordered <- validate_window_phenotype_manifest(manifest, "chr1_0_2000000")
  expect_identical_value(ordered$phenotype_id, c("trait_a", "trait_c", "trait_b"), "phenotype order")
  expect_identical_value(ordered$processing_index, 0:2, "processing indices")
  expect_identical_value(length(unique(ordered$phenotype_key)), 3L, "phenotype key count")
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
    checkpoint_phenotype_key("chr1_0_2000000", "expression", "trait_a"),
    "phenotype key"
  )
})

run_named_test("file hash is deterministic", {
  test_file <- tempfile(fileext = ".txt")
  on.exit(unlink(test_file), add = TRUE)
  writeLines("canonical content", test_file)
  expect_identical_value(sha256_file(test_file), sha256_file(test_file), "file hash")
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
