source("tests/test_helpers.R")

Sys.setenv(
  SUSIER_FUNCTIONS_PATH = normalizePath("R/utils", mustWork = TRUE),
  CHECKPOINTED_SUSIE_BASE_IMAGE_DIGEST = "sha256:controller-integration-test"
)
source("R/scripts/run_checkpointed_window_susie.R")

expect_true_value <- function(value, label) {
  if (!isTRUE(value)) {
    stop("Expected true for ", label, ".", call. = FALSE)
  }
  invisible(TRUE)
}

expect_error_condition <- function(code, class, message = NULL) {
  condition <- tryCatch(
    {
      force(code)
      NULL
    },
    error = identity
  )
  if (is.null(condition)) {
    stop("Expected an error of class ", class, ".", call. = FALSE)
  }
  if (!inherits(condition, class)) {
    stop(
      "Expected error class ", class, " but received ",
      paste(class(condition), collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!is.null(message) && !grepl(message, conditionMessage(condition), fixed = TRUE)) {
    stop(
      "Expected error containing '", message,
      "' but received: ", conditionMessage(condition),
      call. = FALSE
    )
  }
  invisible(condition)
}

fixture_path <- function(name) {
  file.path("tests", "fixtures", "checkpointed_window", name)
}

controller_script <- normalizePath(
  "R/scripts/run_checkpointed_window_susie.R",
  mustWork = TRUE
)

make_controller_config <- function(label) {
  list(
    window_id = "chr1_0_2000000",
    window_dosage = fixture_path("window_dosage.tsv"),
    window_phenotypes = fixture_path("window_phenotypes.tsv"),
    phenotype_data = fixture_path("window_phenotypes.bed.gz"),
    covariate_files = c(
      fixture_path("shared_covariates.tsv"),
      fixture_path("expression_covariates.tsv")
    ),
    covariate_modalities = c("shared", "expression"),
    keep_samples = fixture_path("keep_samples.tsv"),
    checkpoint_root = tempfile(paste0("checkpointed-controller-", label, "-")),
    output_dir = tempfile(paste0("checkpointed-output-", label, "-")),
    wrapper_path = controller_script,
    source_paths = c(
      runner_wrapper = controller_script,
      checkpointed_functions = normalizePath(
        "R/utils/CheckpointedWindowSusieFunctions.R",
        mustWork = TRUE
      ),
      checkpoint_store = normalizePath("R/utils/CheckpointStore.R", mustWork = TRUE)
    ),
    runtime_image = "local-controller-test",
    base_image_digest = Sys.getenv("CHECKPOINTED_SUSIE_BASE_IMAGE_DIGEST")
  )
}

make_three_phenotype_config <- function(label) {
  config <- make_controller_config(label)
  input_dir <- tempfile(paste0("checkpointed-three-phenotype-", label, "-"))
  dir.create(input_dir, recursive = TRUE)
  phenotype_data <- data.table::fread(
    fixture_path("window_phenotypes.bed.gz"),
    check.names = FALSE,
    data.table = FALSE
  ) |>
    tibble::as_tibble()
  third_phenotype <- phenotype_data |>
    dplyr::filter(.data$phenotype_id == "linked_expression") |>
    dplyr::mutate(
      chromosome = 8L,
      start = 800000L,
      end = 800001L,
      phenotype_id = "linked_expression_replica",
      dplyr::across(
        dplyr::starts_with("sample_"),
        ~ .x + seq_along(.x) * 1e-4
      )
    )
  phenotype_path <- file.path(input_dir, "window_phenotypes.tsv")
  readr::write_tsv(
    dplyr::bind_rows(phenotype_data, third_phenotype),
    phenotype_path
  )
  manifest_path <- file.path(input_dir, "manifest.tsv")
  readr::write_tsv(
    tibble::tribble(
      ~window_id, ~phenotype_id, ~modality, ~phenotype_file, ~p_value,
      config$window_id, "linked_expression", "expression", phenotype_path, 1e-12,
      config$window_id, "linked_splicing", "splicing", phenotype_path, 2e-8,
      config$window_id, "linked_expression_replica", "expression", phenotype_path, 3e-8
    ),
    manifest_path
  )
  config$window_phenotypes <- manifest_path
  config$phenotype_data <- phenotype_path
  config$test_input_dir <- input_dir
  config
}

controller_ordered_manifest <- function(config) {
  raw_manifest <- readr::read_tsv(
    config$window_phenotypes,
    show_col_types = FALSE
  )
  ordered_manifest <- validate_window_phenotype_manifest(
    raw_manifest,
    config$window_id
  )
  covariate_hash_names <- paste0(
    "covariate_",
    seq_along(config$covariate_files),
    "_",
    config$covariate_modalities
  )
  scientific_paths <- c(
    window_dosage = config$window_dosage,
    window_phenotypes = config$window_phenotypes,
    phenotype_data = config$phenotype_data,
    stats::setNames(config$covariate_files, covariate_hash_names),
    keep_samples = config$keep_samples
  )
  input_hashes <- purrr::map_chr(scientific_paths, sha256_file)
  analysis_id <- build_checkpoint_analysis_id(
    input_hashes,
    ordered_manifest,
    checkpointed_susie_settings(),
    Sys.getenv("CHECKPOINTED_SUSIE_BASE_IMAGE_DIGEST"),
    sha256_file(config$wrapper_path),
    source_hashes = purrr::map_chr(config$source_paths, sha256_file),
    runtime_image = config$runtime_image
  )
  ordered_manifest |>
    dplyr::mutate(analysis_id = analysis_id, .before = 1L)
}

read_controller_window_manifest <- function(config, store = NULL) {
  ordered <- controller_ordered_manifest(config)
  paths <- window_checkpoint_paths(
    ordered$analysis_id[[1L]],
    config$window_id
  )
  store <- store %||% new_checkpoint_store(config$checkpoint_root)
  read_window_run_manifest_if_present(store, paths$window_manifest)
}

read_controller_fit_manifest <- function(config, processing_index, store = NULL) {
  ordered <- controller_ordered_manifest(config)
  store <- store %||% new_checkpoint_store(config$checkpoint_root)
  read_valid_boundary_manifest(
    store,
    ordered[ordered$processing_index == processing_index, , drop = FALSE]
  )
}

prepared_config <- make_controller_config("prepared-inputs")
prepared_manifest <- controller_ordered_manifest(prepared_config)
prepared_inputs <- load_checkpointed_window_inputs(
  prepared_config,
  prepared_manifest
)
invisible(utils::capture.output(
  prepared_results <- purrr::map(
    seq_len(nrow(prepared_manifest)),
    function(row_index) {
      record <- prepared_manifest[row_index, , drop = FALSE]
      fit_checkpointed_window_phenotype(
        genotype = prepared_inputs$genotype,
        phenotype = prepared_inputs$phenotypes[[record$phenotype_id[[1L]]]],
        covariates = applicable_covariate_matrix(
          prepared_inputs$covariates,
          record$modality[[1L]],
          prepared_inputs$sample_ids
        ),
        variant_ids = prepared_inputs$variant_info$variant_id,
        settings = checkpointed_susie_settings()
      )
    }
  ) |>
    stats::setNames(prepared_manifest$phenotype_id)
))

new_counting_fit_function <- function(
    counts,
    results = prepared_results,
    failure = NULL,
    phenotype_values = prepared_inputs$phenotypes,
    before_fit = NULL
) {
  force(counts)
  force(results)
  force(failure)
  force(phenotype_values)
  force(before_fit)
  function(genotype, phenotype, covariates, variant_ids, settings) {
    phenotype_matches <- purrr::map_lgl(
      phenotype_values,
      ~ identical(as.numeric(.x), as.numeric(phenotype))
    )
    phenotype_id <- names(phenotype_values)[phenotype_matches]
    if (length(phenotype_id) != 1L) {
      stop("The counting fit function could not identify the phenotype.", call. = FALSE)
    }
    if (!is.null(before_fit)) {
      before_fit(phenotype_id)
    }
    counts$by_phenotype[[phenotype_id]] <-
      counts$by_phenotype[[phenotype_id]] + 1L
    if (!is.null(failure)) {
      stop(failure)
    }
    results[[phenotype_id]]
  }
}

new_fit_counts <- function(phenotype_ids = names(prepared_inputs$phenotypes)) {
  counts <- new.env(parent = emptyenv())
  counts$by_phenotype <- stats::setNames(
    rep(0L, length(phenotype_ids)),
    phenotype_ids
  )
  counts
}

prepare_three_phenotype_case <- function(config) {
  ordered <- controller_ordered_manifest(config)
  inputs <- load_checkpointed_window_inputs(config, ordered)
  results <- list(
    linked_expression = prepared_results[["linked_expression"]],
    linked_splicing = prepared_results[["linked_splicing"]],
    linked_expression_replica = prepared_results[["linked_expression"]]
  )
  list(ordered = ordered, inputs = inputs, results = results)
}

commit_direct_fit <- function(config, store, processing_index, fit_result) {
  ordered <- controller_ordered_manifest(config)
  record <- ordered[ordered$processing_index == processing_index, , drop = FALSE]
  fit_result$started_at <- "2026-08-27T12:00:00Z"
  fit_result$completed_at <- "2026-08-27T12:01:00Z"
  local_artifacts <- write_local_phenotype_artifacts(
    fit_result,
    record,
    tempfile("direct-checkpoint-artifacts-")
  )
  covariate_hash_names <- paste0(
    "covariate_",
    seq_along(config$covariate_files),
    "_",
    config$covariate_modalities
  )
  scientific_paths <- c(
    window_dosage = config$window_dosage,
    window_phenotypes = config$window_phenotypes,
    phenotype_data = config$phenotype_data,
    stats::setNames(config$covariate_files, covariate_hash_names),
    keep_samples = config$keep_samples
  )
  input_hashes <- purrr::map_chr(scientific_paths, sha256_file)
  fit_manifest <- build_fit_manifest(
    ordered$analysis_id[[1L]],
    config$window_id,
    record,
    fit_result,
    local_artifacts,
    input_hashes,
    checkpointed_susie_settings(),
    checkpointed_runtime_provenance(config)
  )
  fit_sha256 <- if (identical(fit_result$status, "SKIPPED")) {
    NULL
  } else {
    fit_manifest$payloads$susie_fit$sha256
  }
  commit_phenotype_checkpoint(
    store,
    phenotype_checkpoint_paths(
      ordered$analysis_id[[1L]],
      config$window_id,
      record$phenotype_key[[1L]],
      fit_sha256
    ),
    local_artifacts,
    fit_manifest
  )
}

run_named_test("controller keeps exact subnormal p values in checkpoint identity", {
  config <- make_controller_config("subnormal-p-values")
  config$window_phenotypes <- fixture_path("subnormal_window_phenotypes.tsv")
  on.exit(unlink(c(config$checkpoint_root, config$output_dir), recursive = TRUE), add = TRUE)
  store <- new_checkpoint_store(config$checkpoint_root)
  counts <- new_fit_counts()
  fit_function <- new_counting_fit_function(counts)
  exact_p_values <- as.numeric(c(
    "1.35242422104231e-309",
    "6.834591488289331e-295"
  ))
  expected_ordered <- tibble::tibble(
    window_id = rep(config$window_id, 2L),
    phenotype_id = c("linked_expression", "linked_splicing"),
    modality = c("expression", "splicing"),
    phenotype_file = rep("window_phenotypes.bed.gz", 2L),
    p_value = exact_p_values,
    processing_index = 0:1,
    phenotype_key = c(
      checkpoint_phenotype_key(
        config$window_id,
        "expression",
        "linked_expression"
      ),
      checkpoint_phenotype_key(
        config$window_id,
        "splicing",
        "linked_splicing"
      )
    )
  )
  input_hashes <- purrr::map_chr(
    checkpointed_scientific_paths(config),
    sha256_file
  )
  expected_analysis_id <- build_checkpoint_analysis_id(
    input_hashes,
    expected_ordered,
    checkpointed_susie_settings(),
    Sys.getenv("CHECKPOINTED_SUSIE_BASE_IMAGE_DIGEST"),
    sha256_file(config$wrapper_path),
    source_hashes = purrr::map_chr(config$source_paths, sha256_file),
    runtime_image = config$runtime_image
  )

  paths <- run_checkpointed_window(
    config,
    store = store,
    fit_function = fit_function
  )
  first_manifest <- jsonlite::read_json(
    paths$window_manifest,
    simplifyVector = FALSE
  )
  expect_identical_value(
    purrr::map_chr(first_manifest$phenotypes, "phenotype_id"),
    expected_ordered$phenotype_id,
    "subnormal-p-value processing order"
  )
  expect_identical_value(
    purrr::map_dbl(first_manifest$phenotypes, "p_value"),
    exact_p_values,
    "checkpoint p values"
  )
  expect_identical_value(
    first_manifest$analysis_id,
    expected_analysis_id,
    "exact-p-value analysis ID"
  )

  resumed_paths <- run_checkpointed_window(
    config,
    store = store,
    fit_function = fit_function
  )
  resumed_manifest <- jsonlite::read_json(
    resumed_paths$window_manifest,
    simplifyVector = FALSE
  )
  expect_identical_value(
    resumed_manifest$analysis_id,
    first_manifest$analysis_id,
    "stable resumed analysis ID"
  )
  expect_identical_value(
    counts$by_phenotype,
    c(linked_expression = 1L, linked_splicing = 1L),
    "resume fit counts"
  )
})

run_named_test("window cursor reader rejects a non-object JSON value", {
  store_root <- tempfile("checkpointed-invalid-window-json-")
  on.exit(unlink(store_root, recursive = TRUE), add = TRUE)
  store <- new_checkpoint_store(store_root)
  local_json <- tempfile(fileext = ".json")
  on.exit(unlink(local_json), add = TRUE)
  writeLines('"not a window manifest"', local_json)
  relative_path <- file.path("analysis", "window", "window_manifest.json")
  store$upload(local_json, relative_path)
  expect_identical_value(
    read_window_run_manifest_if_present(store, relative_path),
    NULL,
    "invalid window JSON"
  )
})

run_named_test("interruption keeps a running cursor and resume fits only the suffix", {
  config <- make_controller_config("interrupt-resume")
  on.exit(unlink(c(config$checkpoint_root, config$output_dir), recursive = TRUE), add = TRUE)
  store <- new_checkpoint_store(config$checkpoint_root)
  counts <- new_fit_counts()
  fit_function <- new_counting_fit_function(counts)

  expect_error_condition(
    run_checkpointed_window(
      config,
      store = store,
      fit_function = fit_function,
      interrupt_after_commits = 1L
    ),
    "checkpoint_test_interrupt",
    "Synthetic interruption"
  )
  interrupted <- read_controller_window_manifest(config, store)
  expect_identical_value(interrupted$status, "RUNNING", "interrupted status")
  expect_identical_value(interrupted$last_committed_index, 0L, "interrupted boundary")
  expect_identical_value(length(interrupted$committed), 1L, "interrupted commits")
  expect_identical_value(
    counts$by_phenotype,
    c(linked_expression = 1L, linked_splicing = 0L),
    "counts after interruption"
  )

  final_paths <- run_checkpointed_window(
    config,
    store = store,
    fit_function = fit_function
  )
  expect_identical_value(
    counts$by_phenotype,
    c(linked_expression = 1L, linked_splicing = 1L),
    "counts after resume"
  )
  expect_identical_value(
    names(final_paths),
    c("window_manifest", "fit_index", "credible_sets", "lbf_variable", "full_susie"),
    "final path names"
  )
  expect_true_value(all(file.exists(unlist(final_paths))), "completed local outputs")
})

run_named_test("injected fit functions receive legacy-imputed retained dosage", {
  config <- make_controller_config("injected-imputed-dosage")
  dosage_path <- tempfile("injected-missing-dosage-", fileext = ".tsv")
  dosage <- readr::read_tsv(
    fixture_path("window_dosage.tsv"),
    show_col_types = FALSE
  )
  dosage$sample_05[[1L]] <- -1
  readr::write_tsv(dosage, dosage_path)
  config$window_dosage <- dosage_path
  on.exit(
    unlink(c(config$checkpoint_root, config$output_dir, dosage_path), recursive = TRUE),
    add = TRUE
  )
  ordered <- controller_ordered_manifest(config)
  raw_inputs <- load_checkpointed_window_inputs(config, ordered)
  expect_true_value(anyNA(raw_inputs$genotype), "raw cache dosage contains missing value")

  observed_missing <- logical()
  counting_fit <- new_counting_fit_function(new_fit_counts())
  injected_fit <- function(genotype, phenotype, covariates, variant_ids, settings) {
    observed_missing <<- c(observed_missing, anyNA(genotype))
    counting_fit(genotype, phenotype, covariates, variant_ids, settings)
  }
  run_checkpointed_window(
    config,
    fit_function = injected_fit
  )
  expect_identical_value(
    observed_missing,
    c(FALSE, FALSE),
    "missing dosage seen by injected fits"
  )
})

run_named_test("resume recomputes a corrupt boundary RDS", {
  config <- make_controller_config("corrupt-boundary")
  on.exit(unlink(c(config$checkpoint_root, config$output_dir), recursive = TRUE), add = TRUE)
  store <- new_checkpoint_store(config$checkpoint_root)
  counts <- new_fit_counts()
  fit_function <- new_counting_fit_function(counts)

  expect_error_condition(
    run_checkpointed_window(
      config,
      store = store,
      fit_function = fit_function,
      interrupt_after_commits = 1L
    ),
    "checkpoint_test_interrupt"
  )
  boundary_manifest <- read_controller_fit_manifest(config, 0L, store)
  writeBin(
    charToRaw("corrupt checkpoint"),
    file.path(config$checkpoint_root, boundary_manifest$payloads$susie_fit$path)
  )

  run_checkpointed_window(config, store = store, fit_function = fit_function)
  expect_identical_value(
    counts$by_phenotype,
    c(linked_expression = 2L, linked_splicing = 1L),
    "counts after corrupt-boundary resume"
  )
  completed <- read_controller_window_manifest(config, store)
  expect_identical_value(length(completed$recovery_history), 1L, "recovery audit count")
  expect_identical_value(
    completed$recovery_history[[1L]]$processing_index,
    0L,
    "recovery audit index"
  )
  expect_true_value(
    grepl("checksum", completed$recovery_history[[1L]]$reason, ignore.case = TRUE),
    "recovery audit reason"
  )
})

run_named_test("final assembly recomputes a corrupt interior Parquet and suffix", {
  config <- make_controller_config("corrupt-interior-parquet")
  on.exit(unlink(c(config$checkpoint_root, config$output_dir), recursive = TRUE), add = TRUE)
  store <- new_checkpoint_store(config$checkpoint_root)
  run_checkpointed_window(
    config,
    store = store,
    fit_function = new_counting_fit_function(new_fit_counts())
  )
  first_manifest <- read_controller_fit_manifest(config, 0L, store)
  corrupt_path <- first_manifest$payloads$full_susie$path
  corrupt_bytes <- charToRaw("corrupt interior Parquet")
  writeBin(
    corrupt_bytes,
    file.path(config$checkpoint_root, corrupt_path)
  )

  retry_counts <- new_fit_counts()
  paths <- run_checkpointed_window(
    config,
    store = store,
    fit_function = new_counting_fit_function(retry_counts)
  )
  expect_identical_value(
    retry_counts$by_phenotype,
    c(linked_expression = 1L, linked_splicing = 1L),
    "corrupt interior recomputation counts"
  )
  completed <- jsonlite::read_json(paths$window_manifest, simplifyVector = FALSE)
  expect_identical_value(completed$status, "COMPLETE", "corruption recovery status")
  expect_identical_value(length(completed$recovery_history), 1L, "corruption audit count")
  expect_identical_value(
    completed$recovery_history[[1L]]$processing_index,
    0L,
    "corrupt Parquet audit index"
  )
  expect_true_value(
    grepl(corrupt_path, completed$recovery_history[[1L]]$reason, fixed = TRUE),
    "corrupt Parquet audit path"
  )
  recovery <- completed$recovery_history[[1L]]
  expect_identical_value(recovery$payload_path, corrupt_path, "recovery source path")
  expect_identical_value(recovery$payload_missing, FALSE, "recovery source presence")
  expect_true_value(nzchar(recovery$evidence_path), "recovery evidence path")
  expected_evidence_directory <- file.path(
    first_manifest$analysis_id,
    config$window_id,
    "recovery_evidence",
    "attempt-00000002",
    "index-00000000"
  )
  expect_identical_value(
    dirname(recovery$evidence_path),
    expected_evidence_directory,
    "recovery evidence directory"
  )
  expect_true_value(
    grepl(
      "^full_susie\\.parquet\\.[0-9a-f]{64}\\.invalid$",
      basename(recovery$evidence_path)
    ),
    "recovery evidence basename"
  )
  expect_identical_value(
    recovery$evidence_uri,
    store$object_uri(recovery$evidence_path),
    "recovery evidence URI"
  )
  evidence_file <- file.path(config$checkpoint_root, recovery$evidence_path)
  evidence_bytes <- readBin(
    evidence_file,
    what = "raw",
    n = unname(file.info(evidence_file)$size)
  )
  expect_identical_value(evidence_bytes, corrupt_bytes, "preserved corrupt bytes")
  repaired_manifest <- read_controller_fit_manifest(config, 0L, store)
  expect_identical_value(
    sha256_file(file.path(config$checkpoint_root, corrupt_path)),
    repaired_manifest$payloads$full_susie$sha256,
    "repaired fixed payload"
  )
  expect_true_value(
    all(file.exists(unlist(paths))),
    "outputs after corruption recovery"
  )
})

run_named_test("final hydration recomputes a missing interior saved fit", {
  config <- make_controller_config("missing-interior-rds")
  on.exit(unlink(c(config$checkpoint_root, config$output_dir), recursive = TRUE), add = TRUE)
  store <- new_checkpoint_store(config$checkpoint_root)
  run_checkpointed_window(
    config,
    store = store,
    fit_function = new_counting_fit_function(new_fit_counts())
  )
  first_manifest <- read_controller_fit_manifest(config, 0L, store)
  missing_path <- first_manifest$payloads$susie_fit$path
  unlink(file.path(config$checkpoint_root, missing_path))
  expect_true_value(!store$object_exists(missing_path), "missing saved fit setup")

  retry_counts <- new_fit_counts()
  paths <- run_checkpointed_window(
    config,
    store = store,
    fit_function = new_counting_fit_function(retry_counts)
  )
  expect_identical_value(
    retry_counts$by_phenotype,
    c(linked_expression = 1L, linked_splicing = 1L),
    "missing saved-fit recomputation counts"
  )
  completed <- jsonlite::read_json(paths$window_manifest, simplifyVector = FALSE)
  expect_identical_value(completed$status, "COMPLETE", "missing saved-fit recovery status")
  expect_identical_value(
    completed$recovery_history[[1L]]$processing_index,
    0L,
    "missing saved-fit audit index"
  )
  expect_identical_value(
    completed$recovery_history[[1L]]$payload_path,
    missing_path,
    "missing saved-fit source path"
  )
  expect_identical_value(
    completed$recovery_history[[1L]]$payload_missing,
    TRUE,
    "missing saved-fit evidence state"
  )
  expect_identical_value(
    completed$recovery_history[[1L]]$evidence_path,
    NULL,
    "missing saved-fit evidence path"
  )
  expect_identical_value(
    completed$recovery_history[[1L]]$evidence_uri,
    NULL,
    "missing saved-fit evidence URI"
  )
  expect_true_value(
    !dir.exists(file.path(
      config$checkpoint_root,
      first_manifest$analysis_id,
      config$window_id,
      "recovery_evidence"
    )),
    "missing saved-fit evidence directory"
  )
  expect_true_value(store$object_exists(missing_path), "recomputed saved fit")
})

run_named_test("evidence copy failure stops before payload recomputation", {
  config <- make_controller_config("failed-evidence-copy")
  on.exit(unlink(c(config$checkpoint_root, config$output_dir), recursive = TRUE), add = TRUE)
  store <- new_checkpoint_store(config$checkpoint_root)
  run_checkpointed_window(
    config,
    store = store,
    fit_function = new_counting_fit_function(new_fit_counts())
  )
  first_manifest <- read_controller_fit_manifest(config, 0L, store)
  corrupt_path <- first_manifest$payloads$full_susie$path
  corrupt_bytes <- charToRaw("corrupt payload before copy failure")
  writeBin(corrupt_bytes, file.path(config$checkpoint_root, corrupt_path))

  failing_store <- store
  failing_store$copy_object <- function(source_path, destination_path) {
    stop_checkpoint_store_error("Synthetic recovery evidence copy failure.")
  }
  retry_counts <- new_fit_counts()
  expect_error_condition(
    run_checkpointed_window(
      config,
      store = failing_store,
      fit_function = new_counting_fit_function(retry_counts)
    ),
    "checkpoint_store_error",
    "Synthetic recovery evidence copy failure."
  )
  expect_identical_value(
    retry_counts$by_phenotype,
    c(linked_expression = 0L, linked_splicing = 0L),
    "copy failure recomputation counts"
  )
  fixed_bytes <- readBin(
    file.path(config$checkpoint_root, corrupt_path),
    what = "raw",
    n = length(corrupt_bytes)
  )
  expect_identical_value(fixed_bytes, corrupt_bytes, "copy failure fixed payload")
})

run_named_test("mid-window resume downloads only the latest saved fit RDS", {
  config <- make_three_phenotype_config("resume-rds-download-count")
  on.exit(
    unlink(
      c(config$checkpoint_root, config$output_dir, config$test_input_dir),
      recursive = TRUE
    ),
    add = TRUE
  )
  case <- prepare_three_phenotype_case(config)
  store <- new_checkpoint_store(config$checkpoint_root)
  initial_counts <- new_fit_counts(names(case$inputs$phenotypes))
  expect_error_condition(
    run_checkpointed_window(
      config,
      store = store,
      fit_function = new_counting_fit_function(
        initial_counts,
        results = case$results,
        phenotype_values = case$inputs$phenotypes
      ),
      interrupt_after_commits = 2L
    ),
    "checkpoint_test_interrupt"
  )
  latest_manifest <- read_controller_fit_manifest(config, 1L, store)

  downloads <- new.env(parent = emptyenv())
  downloads$paths <- character()
  recording_store <- store
  recording_store$download <- function(relative_path, local_path) {
    downloads$paths <- c(downloads$paths, relative_path)
    store$download(relative_path, local_path)
  }
  rds_before_fit <- character()
  retry_counts <- new_fit_counts(names(case$inputs$phenotypes))
  run_checkpointed_window(
    config,
    store = recording_store,
    fit_function = new_counting_fit_function(
      retry_counts,
      results = case$results,
      phenotype_values = case$inputs$phenotypes,
      before_fit = function(phenotype_id) {
        if (identical(phenotype_id, "linked_expression_replica")) {
          rds_before_fit <<- downloads$paths[endsWith(downloads$paths, ".rds")]
        }
      }
    )
  )
  expect_identical_value(
    rds_before_fit,
    latest_manifest$payloads$susie_fit$path,
    "latest RDS before resumed fit"
  )
  expect_identical_value(
    downloads$paths[endsWith(downloads$paths, ".rds")],
    latest_manifest$payloads$susie_fit$path,
    "total RDS downloads during resume"
  )
  expect_true_value(
    sum(endsWith(downloads$paths, ".parquet")) >= 9L,
    "final Parquet downloads"
  )
})

run_named_test("phenotype manifests contain cache QC and runtime provenance", {
  config <- make_controller_config("manifest-readiness")
  on.exit(unlink(c(config$checkpoint_root, config$output_dir), recursive = TRUE), add = TRUE)
  store <- new_checkpoint_store(config$checkpoint_root)
  invisible(utils::capture.output(
    run_checkpointed_window(config, store = store)
  ))
  fit_manifest <- read_controller_fit_manifest(config, 0L, store)
  window_manifest <- read_controller_window_manifest(config, store)

  expect_identical_value(fit_manifest$counts$raw_dosage_samples, 40L, "raw samples")
  expect_identical_value(fit_manifest$counts$overlap_samples, 36L, "overlap samples")
  expect_identical_value(
    fit_manifest$counts$phenotype_retained_samples,
    36L,
    "phenotype-retained samples"
  )
  expect_true_value(nzchar(fit_manifest$qc$design_id), "fit design identity")
  expect_true_value(nzchar(fit_manifest$qc$cache_key), "fit cache key")
  expect_identical_value(
    names(fit_manifest$source_hashes),
    names(config$source_paths),
    "fit source hashes"
  )
  expect_identical_value(fit_manifest$runtime$runner_image, config$runtime_image, "fit image")
  expect_identical_value(
    names(window_manifest$source_hashes),
    names(config$source_paths),
    "window source hashes"
  )
  expect_identical_value(
    window_manifest$runtime$base_image_digest,
    config$base_image_digest,
    "window base image"
  )
})

run_named_test("controller logs pre-fit cache identity and sample count", {
  config <- make_controller_config("pre-fit-log")
  on.exit(unlink(c(config$checkpoint_root, config$output_dir), recursive = TRUE), add = TRUE)
  messages <- capture.output(
    invisible(run_checkpointed_window(
      config,
      fit_function = new_counting_fit_function(new_fit_counts())
    )),
    type = "message"
  )
  expect_true_value(
    any(grepl(
      "Fit index 0 phenotype linked_expression modality expression samples 36 design",
      messages,
      fixed = TRUE
    )),
    "pre-fit log fields"
  )
})

run_named_test("lagging adoption uploads RUNNING before the following fit", {
  config <- make_three_phenotype_config("lagging-cursor")
  on.exit(
    unlink(
      c(config$checkpoint_root, config$output_dir, config$test_input_dir),
      recursive = TRUE
    ),
    add = TRUE
  )
  case <- prepare_three_phenotype_case(config)
  expect_identical_value(nrow(case$ordered), 3L, "lagging test phenotype count")
  store <- new_checkpoint_store(config$checkpoint_root)
  counts <- new_fit_counts(names(case$inputs$phenotypes))
  fit_function <- new_counting_fit_function(
    counts,
    results = case$results,
    phenotype_values = case$inputs$phenotypes
  )

  expect_error_condition(
    run_checkpointed_window(
      config,
      store = store,
      fit_function = fit_function,
      interrupt_after_commits = 1L
    ),
    "checkpoint_test_interrupt"
  )
  commit_direct_fit(
    config,
    store,
    1L,
    case$results[["linked_splicing"]]
  )

  before_fit <- function(phenotype_id) {
    expect_identical_value(
      phenotype_id,
      "linked_expression_replica",
      "post-adoption phenotype"
    )
    adopted_cursor <- read_controller_window_manifest(config, store)
    expect_identical_value(adopted_cursor$status, "RUNNING", "adopted cursor status")
    expect_identical_value(adopted_cursor$last_committed_index, 1L, "adopted cursor boundary")
  }
  run_checkpointed_window(
    config,
    store = store,
    fit_function = new_counting_fit_function(
      counts,
      results = case$results,
      phenotype_values = case$inputs$phenotypes,
      before_fit = before_fit
    )
  )
  expect_identical_value(
    counts$by_phenotype,
    c(
      linked_expression = 1L,
      linked_splicing = 0L,
      linked_expression_replica = 1L
    ),
    "counts after lagging-cursor adoption"
  )
  completed <- read_controller_window_manifest(config, store)
  expect_identical_value(completed$status, "COMPLETE", "adopted completion status")
  expect_identical_value(completed$last_committed_index, 2L, "adopted final boundary")
})

run_named_test("the fit index has exact values for every terminal state", {
  config <- make_three_phenotype_config("fit-index-schema")
  on.exit(
    unlink(
      c(config$checkpoint_root, config$output_dir, config$test_input_dir),
      recursive = TRUE
    ),
    add = TRUE
  )
  case <- prepare_three_phenotype_case(config)
  terminal_results <- case$results
  terminal_results[["linked_expression"]]$status <- "CONVERGED"
  terminal_results[["linked_expression"]]$fit$converged <- TRUE
  terminal_results[["linked_splicing"]] <- list(
    status = "SKIPPED",
    fit = NULL,
    reason = "NO_USABLE_VARIANTS",
    qc = list(
      input_samples = 36L,
      retained_samples = 36L,
      input_variants = 6L,
      retained_variants = 0L,
      removed_variant_ids = case$inputs$variant_info$variant_id
    )
  )
  terminal_results[["linked_expression_replica"]]$status <- "NONCONVERGED"
  terminal_results[["linked_expression_replica"]]$fit$converged <- FALSE
  store <- new_checkpoint_store(config$checkpoint_root)

  paths <- run_checkpointed_window(
    config,
    store = store,
    fit_function = new_counting_fit_function(
      new_fit_counts(names(case$inputs$phenotypes)),
      results = terminal_results,
      phenotype_values = case$inputs$phenotypes
    )
  )
  fit_manifests <- purrr::map(
    0:2,
    ~ read_controller_fit_manifest(config, .x, store)
  )
  fitted_paths <- purrr::map_chr(
    fit_manifests[c(1L, 3L)],
    ~ .x$payloads$susie_fit$path
  )
  fitted_sha256 <- purrr::map_chr(
    file.path(config$checkpoint_root, fitted_paths),
    sha256_file
  )
  fit_index <- readr::read_tsv(paths$fit_index, show_col_types = FALSE) |>
    dplyr::mutate(processing_index = as.integer(.data$processing_index))
  expected_fit_index <- tibble::tibble(
    window_id = rep(config$window_id, 3L),
    processing_index = 0:2,
    phenotype_id = case$ordered$phenotype_id,
    phenotype_key = case$ordered$phenotype_key,
    modality = case$ordered$modality,
    p_value = case$ordered$p_value,
    status = c("CONVERGED", "SKIPPED", "NONCONVERGED"),
    converged = c(TRUE, NA, FALSE),
    exclusion_reason = c(NA_character_, "NO_USABLE_VARIANTS", NA_character_),
    fit_manifest_uri = purrr::map_chr(
      fit_manifests,
      ~ store$object_uri(.x$fit_manifest_path)
    ),
    susie_fit_uri = c(
      store$object_uri(fitted_paths[[1L]]),
      NA_character_,
      store$object_uri(fitted_paths[[2L]])
    ),
    susie_fit_sha256 = c(
      fitted_sha256[[1L]],
      NA_character_,
      fitted_sha256[[2L]]
    )
  )

  expect_identical_value(
    names(fit_index),
    names(expected_fit_index),
    "fit-index columns"
  )
  expect_identical_value(fit_index, expected_fit_index, "complete fit-index schema")
  completed <- read_controller_window_manifest(config, store)
  expect_identical_value(completed$last_committed_index, 2L, "terminal-state boundary")
  expect_identical_value(
    readRDS(file.path(config$checkpoint_root, fitted_paths[[2L]]))$converged,
    FALSE,
    "saved nonconverged RDS"
  )
})

run_named_test("an unexpected fit failure writes a failed cursor and rethrows", {
  config <- make_controller_config("unexpected-failure")
  on.exit(unlink(c(config$checkpoint_root, config$output_dir), recursive = TRUE), add = TRUE)
  counts <- new_fit_counts()
  store <- new_checkpoint_store(config$checkpoint_root)

  expect_error_condition(
    run_checkpointed_window(
      config,
      store = store,
      fit_function = new_counting_fit_function(
        counts,
        failure = simpleError("synthetic unexpected failure")
      )
    ),
    "simpleError",
    "synthetic unexpected failure"
  )
  failed <- read_controller_window_manifest(config, store)
  expect_identical_value(failed$status, "FAILED", "unexpected failure status")
  expect_identical_value(failed$failure$processing_index, 0L, "failure index")
  expect_identical_value(
    failed$failure$phenotype_id,
    "linked_expression",
    "failure phenotype"
  )
  expect_identical_value(
    failed$failure$message,
    "synthetic unexpected failure",
    "failure message"
  )
})

run_named_test("resume preserves prior failure history and clears the active failure", {
  config <- make_controller_config("failure-history")
  on.exit(unlink(c(config$checkpoint_root, config$output_dir), recursive = TRUE), add = TRUE)
  store <- new_checkpoint_store(config$checkpoint_root)
  failing_counts <- new_fit_counts()
  expect_error_condition(
    run_checkpointed_window(
      config,
      store = store,
      fit_function = new_counting_fit_function(
        failing_counts,
        failure = simpleError("prior failed attempt")
      )
    ),
    "simpleError"
  )

  retry_counts <- new_fit_counts()
  paths <- run_checkpointed_window(
    config,
    store = store,
    fit_function = new_counting_fit_function(retry_counts)
  )
  completed <- jsonlite::read_json(paths$window_manifest, simplifyVector = FALSE)
  expect_identical_value(completed$failure, NULL, "active failure after resume")
  expect_identical_value(length(completed$failure_history), 1L, "failure history length")
  expect_identical_value(
    completed$failure_history[[1L]]$message,
    "prior failed attempt",
    "preserved failure message"
  )
})

run_named_test("fallback recovery hydrates an exact-path prefix for final assembly", {
  config <- make_controller_config("recovered-prefix")
  on.exit(unlink(c(config$checkpoint_root, config$output_dir), recursive = TRUE), add = TRUE)
  store <- new_checkpoint_store(config$checkpoint_root)
  commit_direct_fit(
    config,
    store,
    0L,
    prepared_results[["linked_expression"]]
  )
  counts <- new_fit_counts()

  paths <- run_checkpointed_window(
    config,
    store = store,
    fit_function = new_counting_fit_function(counts)
  )
  expect_identical_value(
    counts$by_phenotype,
    c(linked_expression = 0L, linked_splicing = 1L),
    "counts after prefix recovery"
  )
  fit_index <- readr::read_tsv(paths$fit_index, show_col_types = FALSE)
  expect_identical_value(
    fit_index$phenotype_id,
    c("linked_expression", "linked_splicing"),
    "recovered fit-index phenotype order"
  )
  completed <- jsonlite::read_json(paths$window_manifest, simplifyVector = FALSE)
  expect_identical_value(length(completed$committed), 2L, "hydrated committed count")
  expect_identical_value(
    purrr::map_int(completed$committed, "processing_index"),
    0:1,
    "hydrated committed order"
  )
})

run_named_test("final hydration rejects changed provenance in an interior manifest", {
  config <- make_controller_config("interior-provenance")
  on.exit(unlink(c(config$checkpoint_root, config$output_dir), recursive = TRUE), add = TRUE)
  store <- new_checkpoint_store(config$checkpoint_root)
  ordered <- controller_ordered_manifest(config)
  first_manifest_path <- checkpoint_fixed_manifest_path(ordered[1L, , drop = FALSE])
  counts <- new_fit_counts()
  fit_function <- new_counting_fit_function(
    counts,
    before_fit = function(phenotype_id) {
      if (!identical(phenotype_id, "linked_splicing")) {
        return(invisible(NULL))
      }
      local_manifest <- file.path(config$checkpoint_root, first_manifest_path)
      manifest <- jsonlite::read_json(local_manifest, simplifyVector = FALSE)
      manifest$source_hashes$checkpoint_store <- paste(rep("9", 64L), collapse = "")
      jsonlite::write_json(
        manifest,
        local_manifest,
        auto_unbox = TRUE,
        pretty = TRUE,
        null = "null",
        digits = NA
      )
    }
  )

  expect_error_message(
    run_checkpointed_window(
      config,
      store = store,
      fit_function = fit_function
    ),
    "Checkpoint source hashes do not match the current analysis."
  )
})

run_named_test("recovered-prefix cursor survives failure and a second interruption", {
  config <- make_three_phenotype_config("recovered-prefix-preemption")
  on.exit(
    unlink(
      c(config$checkpoint_root, config$output_dir, config$test_input_dir),
      recursive = TRUE
    ),
    add = TRUE
  )
  case <- prepare_three_phenotype_case(config)
  store <- new_checkpoint_store(config$checkpoint_root)
  commit_direct_fit(
    config,
    store,
    0L,
    case$results[["linked_expression"]]
  )

  failed_counts <- new_fit_counts(names(case$inputs$phenotypes))
  expect_error_condition(
    run_checkpointed_window(
      config,
      store = store,
      fit_function = new_counting_fit_function(
        failed_counts,
        results = case$results,
        failure = simpleError("recovered-prefix failed attempt"),
        phenotype_values = case$inputs$phenotypes
      )
    ),
    "simpleError",
    "recovered-prefix failed attempt"
  )
  failed <- read_controller_window_manifest(config, store)
  expect_identical_value(failed$status, "FAILED", "recovered-prefix failed status")
  expect_identical_value(
    failed$recovered_prefix_last_committed_index,
    0L,
    "recovered-prefix boundary"
  )
  expect_identical_value(length(failed$committed), 0L, "recovered-prefix empty suffix")

  retry_counts <- new_fit_counts(names(case$inputs$phenotypes))
  retry_fit <- new_counting_fit_function(
    retry_counts,
    results = case$results,
    phenotype_values = case$inputs$phenotypes
  )
  expect_error_condition(
    run_checkpointed_window(
      config,
      store = store,
      fit_function = retry_fit,
      interrupt_after_commits = 1L
    ),
    "checkpoint_test_interrupt"
  )
  interrupted <- read_controller_window_manifest(config, store)
  expect_identical_value(interrupted$status, "RUNNING", "second interruption status")
  expect_identical_value(interrupted$last_committed_index, 1L, "second interruption boundary")
  expect_identical_value(
    purrr::map_int(interrupted$committed, "processing_index"),
    1L,
    "recovered-prefix suffix inventory"
  )
  expect_identical_value(length(interrupted$failure_history), 1L, "interrupted failure history")

  resumable <- resolve_resume_boundary(store, case$ordered, interrupted)
  expect_identical_value(resumable$next_index, 2L, "partial cursor next index")
  expect_true_value(!is.null(resumable$window_manifest), "partial cursor reuse")

  paths <- run_checkpointed_window(
    config,
    store = store,
    fit_function = retry_fit
  )
  expect_identical_value(
    retry_counts$by_phenotype,
    c(
      linked_expression = 0L,
      linked_splicing = 1L,
      linked_expression_replica = 1L
    ),
    "counts across recovered-prefix interruptions"
  )
  completed <- jsonlite::read_json(paths$window_manifest, simplifyVector = FALSE)
  expect_identical_value(completed$status, "COMPLETE", "recovered-prefix final status")
  expect_identical_value(length(completed$failure_history), 1L, "final failure history")
  expect_identical_value(
    purrr::map_int(completed$committed, "processing_index"),
    0:2,
    "hydrated complete inventory"
  )
})

run_named_test("payloads precede fit manifests and cursors in durable upload order", {
  config <- make_controller_config("upload-order")
  on.exit(unlink(c(config$checkpoint_root, config$output_dir), recursive = TRUE), add = TRUE)
  filesystem_store <- new_checkpoint_store(config$checkpoint_root)
  uploads <- new.env(parent = emptyenv())
  uploads$paths <- character()
  recording_store <- filesystem_store
  recording_store$upload <- function(local_path, relative_path) {
    uploads$paths <- c(uploads$paths, relative_path)
    filesystem_store$upload(local_path, relative_path)
  }

  run_checkpointed_window(
    config,
    store = recording_store,
    fit_function = new_counting_fit_function(new_fit_counts())
  )
  ordered <- controller_ordered_manifest(config)
  window_paths <- window_checkpoint_paths(
    ordered$analysis_id[[1L]],
    config$window_id
  )
  phenotype_uploads <- purrr::map(
    0:1,
    function(processing_index) {
      manifest <- read_controller_fit_manifest(
        config,
        processing_index,
        filesystem_store
      )
      c(
        purrr::map_chr(
          manifest$payloads[c(
            "susie_fit", "credible_sets", "lbf_variable", "full_susie"
          )],
          "path"
        ),
        manifest$fit_manifest_path,
        window_paths$window_manifest
      )
    }
  ) |>
    unlist(use.names = FALSE)
  expected_uploads <- c(
    phenotype_uploads,
    unname(unlist(window_paths[c(
      "fit_index", "credible_sets", "lbf_variable", "full_susie"
    )])),
    window_paths$window_manifest
  )
  expect_identical_value(uploads$paths, expected_uploads, "durable upload order")
})

run_named_test("CLI accepts all flags and rejects invalid covariate arrays", {
  args <- c(
    "--window-id", "window",
    "--window-dosage", "dosage.tsv",
    "--window-phenotypes", "manifest.tsv",
    "--phenotype-data", "phenotypes.bed.gz",
    "--covariate-files", "shared.tsv,expression.tsv",
    "--covariate-modalities", "shared,expression",
    "--keep-samples", "keep.tsv",
    "--checkpoint-root", "gs://bucket/checkpoints",
    "--output-dir", "output"
  )
  parsed <- parse_checkpointed_window_cli(args)
  expect_identical_value(
    names(parsed),
    c(
      "window_id", "window_dosage", "window_phenotypes", "phenotype_data",
      "covariate_files", "covariate_modalities", "keep_samples",
      "checkpoint_root", "output_dir", "wrapper_path", "source_paths",
      "runtime_image", "base_image_digest"
    ),
    "CLI config fields"
  )
  expect_identical_value(
    parsed$covariate_files,
    c("shared.tsv", "expression.tsv"),
    "CLI covariate files"
  )
  expect_identical_value(
    parsed$covariate_modalities,
    c("shared", "expression"),
    "CLI covariate modalities"
  )

  mismatched <- args
  mismatched[match("shared,expression", mismatched)] <- "shared"
  expect_error_message(
    parse_checkpointed_window_cli(mismatched),
    "equal nonzero lengths"
  )
  empty_member <- args
  empty_member[match("shared.tsv,expression.tsv", empty_member)] <-
    "shared.tsv,,expression.tsv"
  expect_error_message(
    parse_checkpointed_window_cli(empty_member),
    "cannot contain empty members"
  )
})

run_named_test("CLI help exits successfully and shows all nine flags", {
  help_output <- system2(
    "Rscript",
    c(shQuote(controller_script), "--help"),
    stdout = TRUE,
    stderr = TRUE
  )
  help_status <- attr(help_output, "status") %||% 0L
  expect_identical_value(help_status, 0L, "CLI help status")
  expected_flags <- c(
    "--window-id", "--window-dosage", "--window-phenotypes",
    "--phenotype-data", "--covariate-files", "--covariate-modalities",
    "--keep-samples", "--checkpoint-root", "--output-dir"
  )
  expect_true_value(
    all(purrr::map_lgl(expected_flags, ~ any(grepl(.x, help_output, fixed = TRUE)))),
    "CLI help flags"
  )
})
