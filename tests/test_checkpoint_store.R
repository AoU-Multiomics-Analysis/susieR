source("tests/test_helpers.R")
source("R/utils/InitFunctions.R")
source("R/utils/CheckpointedWindowSusieFunctions.R")
source("R/utils/CheckpointStore.R")

expect_true_value <- function(value, label) {
  if (!isTRUE(value)) {
    stop("Expected true for ", label, ".", call. = FALSE)
  }
  invisible(TRUE)
}

fake_fit <- structure(
  list(
    pip = c(0.1, 0.9),
    alpha = matrix(c(0.1, 0.9), nrow = 1L),
    mu = matrix(c(0.2, 0.4), nrow = 1L),
    mu2 = matrix(c(0.05, 0.2), nrow = 1L),
    variant_id = c("chr1_100_A_G", "chr1_200_C_T"),
    converged = TRUE
  ),
  class = "susie"
)

make_ordered_manifest <- function(count, analysis_id = "analysis-1", window_id = "chr1_0_2000000") {
  tibble::tibble(
    analysis_id = rep(analysis_id, count),
    window_id = rep(window_id, count),
    phenotype_id = sprintf("trait_%04d", seq_len(count)),
    modality = rep("expression", count),
    phenotype_file = rep("window.bed.gz", count),
    p_value = seq_len(count) / (count + 1),
    processing_index = seq_len(count) - 1L,
    phenotype_key = sprintf("key_%04d", seq_len(count))
  )
}

make_local_artifacts <- function(fit = fake_fit) {
  artifact_directory <- tempfile("checkpoint-artifacts-")
  dir.create(artifact_directory, recursive = TRUE)
  paths <- list(
    fit_rds = file.path(artifact_directory, "fit.rds"),
    credible_sets = file.path(artifact_directory, "credible_sets.parquet"),
    lbf_variable = file.path(artifact_directory, "lbf_variable.parquet"),
    full_susie = file.path(artifact_directory, "full_susie.parquet")
  )
  saveRDS(fit, paths$fit_rds)
  write_checkpointed_susie_tables(
    checkpointed_empty_susie_tables(),
    artifact_directory
  )
  paths
}

payload_record <- function(local_path, relative_path) {
  list(
    path = relative_path,
    sha256 = sha256_file(local_path),
    bytes = unname(file.info(local_path)$size)
  )
}

make_fit_manifest <- function(record, local_artifacts, paths, fit = fake_fit) {
  list(
    schema_version = 1L,
    analysis_id = record$analysis_id[[1L]],
    window_id = record$window_id[[1L]],
    phenotype_id = record$phenotype_id[[1L]],
    phenotype_key = record$phenotype_key[[1L]],
    modality = record$modality[[1L]],
    processing_index = record$processing_index[[1L]],
    p_value = record$p_value[[1L]],
    status = if (isTRUE(fit$converged)) "CONVERGED" else "NONCONVERGED",
    converged = fit$converged,
    exclusion_reason = NULL,
    input_hashes = list(
      dosage = paste(rep("a", 64L), collapse = ""),
      phenotypes = paste(rep("b", 64L), collapse = "")
    ),
    settings = checkpointed_susie_settings(),
    container_digest = "sha256:checkpointed-window-susie",
    package_versions = list(R = "4.5.1", susieR = "0.12.35"),
    source_hashes = list(
      runner_wrapper = paste(rep("c", 64L), collapse = ""),
      checkpointed_functions = paste(rep("d", 64L), collapse = ""),
      checkpoint_store = paste(rep("e", 64L), collapse = "")
    ),
    runtime = list(
      runner_image = "ghcr.io/example/checkpointed@sha256:published",
      base_image_digest = "sha256:checkpointed-window-susie"
    ),
    counts = list(
      input_samples = 40L,
      retained_samples = 38L,
      raw_dosage_samples = 40L,
      overlap_samples = 40L,
      phenotype_retained_samples = 38L,
      input_variants = length(fit$variant_id),
      retained_variants = length(fit$variant_id)
    ),
    qc = list(
      removed_covariates = character(),
      removed_variant_ids = character(),
      design_id = paste(rep("f", 64L), collapse = ""),
      sample_mask_id = paste(rep("1", 64L), collapse = ""),
      cache_key = paste(rep("0", 64L), collapse = "")
    ),
    started_at = "2026-08-27T12:00:00Z",
    completed_at = "2026-08-27T12:01:00Z",
    payloads = list(
      susie_fit = payload_record(local_artifacts$fit_rds, paths$fit_rds),
      credible_sets = payload_record(local_artifacts$credible_sets, paths$credible_sets),
      lbf_variable = payload_record(local_artifacts$lbf_variable, paths$lbf_variable),
      full_susie = payload_record(local_artifacts$full_susie, paths$full_susie)
    )
  )
}

make_skipped_commit_inputs <- function(record, reason = "NO_USABLE_VARIANTS") {
  local_artifacts <- make_local_artifacts()
  fitted_paths <- phenotype_checkpoint_paths(
    record$analysis_id[[1L]],
    record$window_id[[1L]],
    record$phenotype_key[[1L]],
    sha256_file(local_artifacts$fit_rds)
  )
  fit_manifest <- make_fit_manifest(record, local_artifacts, fitted_paths)
  local_artifacts$fit_rds <- NULL
  paths <- phenotype_checkpoint_paths(
    record$analysis_id[[1L]],
    record$window_id[[1L]],
    record$phenotype_key[[1L]]
  )
  fit_manifest$status <- "SKIPPED"
  fit_manifest["converged"] <- list(NULL)
  fit_manifest$exclusion_reason <- reason
  fit_manifest$payloads <- list(
    credible_sets = payload_record(local_artifacts$credible_sets, paths$credible_sets),
    lbf_variable = payload_record(local_artifacts$lbf_variable, paths$lbf_variable),
    full_susie = payload_record(local_artifacts$full_susie, paths$full_susie)
  )
  list(
    local_artifacts = local_artifacts,
    paths = paths,
    fit_manifest = fit_manifest
  )
}

make_commit_inputs <- function(record, fit = fake_fit) {
  local_artifacts <- make_local_artifacts(fit)
  fit_sha256 <- sha256_file(local_artifacts$fit_rds)
  paths <- phenotype_checkpoint_paths(
    record$analysis_id[[1L]],
    record$window_id[[1L]],
    record$phenotype_key[[1L]],
    fit_sha256
  )
  list(
    local_artifacts = local_artifacts,
    paths = paths,
    fit_manifest = make_fit_manifest(record, local_artifacts, paths, fit)
  )
}

run_named_test("local checkpoint store copies exact paths", {
  store_root <- tempfile("checkpoint-store-")
  store <- new_checkpoint_store(paste0(store_root, "/"))
  local_source <- tempfile(fileext = ".txt")
  local_download <- tempfile(fileext = ".txt")
  writeLines("durable payload", local_source)

  store$upload(local_source, file.path("analysis", "window", "payload.txt"))
  expect_true_value(
    store$object_exists(file.path("analysis", "window", "payload.txt")),
    "uploaded local object"
  )
  store$download(file.path("analysis", "window", "payload.txt"), local_download)
  expect_identical_value(readLines(local_download), "durable payload", "downloaded payload")
  expect_identical_value(store$root, store_root, "clean store root")
})

run_named_test("checkpoint object copy preserves bytes and reports local failure", {
  store_root <- tempfile("checkpoint-copy-")
  store <- new_checkpoint_store(store_root)
  source_path <- file.path("analysis", "window", "payload.bin")
  evidence_path <- file.path("analysis", "window", "recovery", "payload.bin")
  source_bytes <- as.raw(c(0L, 1L, 2L, 127L, 255L))
  local_source <- tempfile(fileext = ".bin")
  writeBin(source_bytes, local_source)
  store$upload(local_source, source_path)

  store$copy_object(source_path, evidence_path)
  copied_bytes <- readBin(
    file.path(store_root, evidence_path),
    what = "raw",
    n = length(source_bytes)
  )
  expect_identical_value(copied_bytes, source_bytes, "copied object bytes")

  condition <- tryCatch(
    {
      store$copy_object("analysis/window/missing.bin", evidence_path)
      NULL
    },
    error = identity
  )
  expect_true_value(
    inherits(condition, "checkpoint_store_error"),
    "local copy error class"
  )
})

run_named_test("checkpoint paths follow the fixed layout", {
  expect_identical_value(
    window_checkpoint_paths("analysis", "window"),
    list(
      window_manifest = file.path("analysis", "window", "window_manifest.json"),
      fit_index = file.path("analysis", "window", "window_fit_index.tsv"),
      credible_sets = file.path("analysis", "window", "results", "credible_sets.parquet"),
      lbf_variable = file.path("analysis", "window", "results", "lbf_variable.parquet"),
      full_susie = file.path("analysis", "window", "results", "full_susie.parquet")
    ),
    "window checkpoint paths"
  )
  expect_identical_value(
    phenotype_checkpoint_paths("analysis", "window", "key", "abc"),
    list(
      fit_manifest = file.path("analysis", "window", "phenotypes", "key", "fit_manifest.json"),
      fit_rds = file.path("analysis", "window", "phenotypes", "key", "susie_fit.abc.rds"),
      credible_sets = file.path("analysis", "window", "phenotypes", "key", "credible_sets.parquet"),
      lbf_variable = file.path("analysis", "window", "phenotypes", "key", "lbf_variable.parquet"),
      full_susie = file.path("analysis", "window", "phenotypes", "key", "full_susie.parquet")
    ),
    "phenotype checkpoint paths"
  )
  expect_identical_value(
    phenotype_checkpoint_paths("analysis", "window", "key")$fit_rds,
    character(0),
    "skipped fit path"
  )
})

run_named_test("fit validator accepts a structurally valid SuSiE fit", {
  expect_true_value(validate_checkpointed_susie_fit(fake_fit), "valid SuSiE fit")
})

run_named_test("fit validator rejects duplicate variant IDs", {
  invalid_fit <- fake_fit
  invalid_fit$variant_id <- rep("chr1_100_A_G", 2L)
  expect_error_message(validate_checkpointed_susie_fit(invalid_fit), "unique")
})

run_named_test("fit validator rejects a PIP outside zero through one", {
  invalid_fit <- fake_fit
  invalid_fit$pip[[1L]] <- 1.1
  expect_error_message(validate_checkpointed_susie_fit(invalid_fit), "between zero and one")
})

run_named_test("fit validator rejects inconsistent matrix dimensions", {
  invalid_fit <- fake_fit
  invalid_fit$mu2 <- matrix(c(0.05, 0.2, 0.3, 0.4), nrow = 2L)
  expect_error_message(validate_checkpointed_susie_fit(invalid_fit), "same dimensions")
})

run_named_test("fit validator rejects nonfinite model arrays", {
  invalid_fit <- fake_fit
  invalid_fit$alpha[[1L]] <- NA_real_
  expect_error_message(validate_checkpointed_susie_fit(invalid_fit), "finite numeric values")
})

run_named_test("fit validator requires vector PIPs and positive model dimensions", {
  matrix_pip <- fake_fit
  matrix_pip$pip <- matrix(matrix_pip$pip, nrow = 1L)
  expect_error_message(validate_checkpointed_susie_fit(matrix_pip), "numeric vector")

  empty_model <- fake_fit
  empty_model$variant_id <- character()
  empty_model$pip <- numeric()
  empty_model$alpha <- matrix(numeric(), nrow = 0L, ncol = 0L)
  empty_model$mu <- matrix(numeric(), nrow = 0L, ncol = 0L)
  empty_model$mu2 <- matrix(numeric(), nrow = 0L, ncol = 0L)
  expect_error_message(validate_checkpointed_susie_fit(empty_model), "positive dimensions")
})

run_named_test("GCS object existence distinguishes not found from operational failure", {
  not_found_gsutil <- tempfile("fake-gsutil-not-found-")
  failed_gsutil <- tempfile("fake-gsutil-failed-")
  on.exit(unlink(c(not_found_gsutil, failed_gsutil)), add = TRUE)
  writeLines(
    c("#!/bin/sh", "echo 'No URLs matched: gs://test-bucket/missing' >&2", "exit 1"),
    not_found_gsutil
  )
  writeLines(
    c("#!/bin/sh", "echo 'AccessDeniedException: quota or authentication failed' >&2", "exit 1"),
    failed_gsutil
  )
  Sys.chmod(c(not_found_gsutil, failed_gsutil), mode = "0755")
  missing_store <- new_checkpoint_store(
    "gs://test-bucket/checkpoints",
    gsutil = not_found_gsutil
  )
  expect_identical_value(
    missing_store$object_exists("analysis/window/missing.json"),
    FALSE,
    "GCS not-found result"
  )

  failed_store <- new_checkpoint_store(
    "gs://test-bucket/checkpoints",
    gsutil = failed_gsutil
  )
  expected_uri <- paste0(
    "gs://test-bucket/checkpoints/analysis/window/manifest.json"
  )
  condition <- tryCatch(
    {
      failed_store$object_exists("analysis/window/manifest.json")
      NULL
    },
    error = identity
  )
  expect_true_value(inherits(condition, "checkpoint_store_error"), "GCS stat error class")
  expect_true_value(
    grepl(expected_uri, conditionMessage(condition), fixed = TRUE),
    "GCS stat error URI"
  )
})

run_named_test("GCS object copy uses server-side URIs and reports operational failure", {
  copy_log <- tempfile("fake-gsutil-copy-log-")
  successful_gsutil <- tempfile("fake-gsutil-copy-")
  failed_gsutil <- tempfile("fake-gsutil-copy-failed-")
  on.exit(unlink(c(copy_log, successful_gsutil, failed_gsutil)), add = TRUE)
  writeLines(
    c(
      "#!/bin/sh",
      paste0("printf '%s\\n' \"$@\" > ", shQuote(copy_log)),
      "exit 0"
    ),
    successful_gsutil
  )
  writeLines(c("#!/bin/sh", "exit 17"), failed_gsutil)
  Sys.chmod(c(successful_gsutil, failed_gsutil), mode = "0755")
  source_path <- "analysis/window/full_susie.parquet"
  evidence_path <- "analysis/window/recovery/full_susie.parquet.invalid"
  source_uri <- paste0("gs://test-bucket/checkpoints/", source_path)
  evidence_uri <- paste0("gs://test-bucket/checkpoints/", evidence_path)

  successful_store <- new_checkpoint_store(
    "gs://test-bucket/checkpoints",
    gsutil = successful_gsutil
  )
  successful_store$copy_object(source_path, evidence_path)
  expect_identical_value(
    readLines(copy_log),
    c("-q", "cp", source_uri, evidence_uri),
    "GCS server-side copy arguments"
  )

  failed_store <- new_checkpoint_store(
    "gs://test-bucket/checkpoints",
    gsutil = failed_gsutil
  )
  condition <- tryCatch(
    {
      failed_store$copy_object(source_path, evidence_path)
      NULL
    },
    error = identity
  )
  expect_true_value(
    inherits(condition, "checkpoint_store_error"),
    "GCS copy error class"
  )
  expect_true_value(
    grepl(source_uri, conditionMessage(condition), fixed = TRUE) &&
      grepl(evidence_uri, conditionMessage(condition), fixed = TRUE),
    "GCS copy error URIs"
  )
})

run_named_test("phenotype commit uploads all payloads before its manifest", {
  record <- make_ordered_manifest(1L)[1L, ]
  commit_inputs <- make_commit_inputs(record)
  calls <- new.env(parent = emptyenv())
  calls$relative_paths <- character()
  mock_store <- list(
    upload = function(local_path, relative_path) {
      if (!file.exists(local_path)) {
        stop("Upload source is absent.", call. = FALSE)
      }
      calls$relative_paths <- c(calls$relative_paths, relative_path)
      invisible(TRUE)
    },
    object_uri = function(relative_path) paste0("gs://test-root/", relative_path)
  )

  committed <- commit_phenotype_checkpoint(
    mock_store,
    commit_inputs$paths,
    commit_inputs$local_artifacts,
    commit_inputs$fit_manifest
  )

  expect_identical_value(
    calls$relative_paths,
    unname(unlist(commit_inputs$paths[c(
      "fit_rds",
      "credible_sets",
      "lbf_variable",
      "full_susie",
      "fit_manifest"
    )])),
    "payload-first upload order"
  )
  expect_identical_value(
    committed$fit_manifest_path,
    commit_inputs$paths$fit_manifest,
    "committed manifest path"
  )
})

run_named_test("checkpoint JSON preserves very small P values for resume", {
  ordered <- make_ordered_manifest(1L) |>
    dplyr::mutate(p_value = 1.1199976554521247e-146)
  store_root <- tempfile("checkpoint-small-p-value-")
  on.exit(unlink(store_root, recursive = TRUE), add = TRUE)
  store <- new_checkpoint_store(store_root)
  commit_inputs <- make_commit_inputs(ordered[1L, ])
  committed <- commit_phenotype_checkpoint(
    store,
    commit_inputs$paths,
    commit_inputs$local_artifacts,
    commit_inputs$fit_manifest
  )

  hydrated_fit_manifest <- read_valid_boundary_manifest(store, ordered[1L, ])
  expect_identical_value(
    hydrated_fit_manifest$p_value,
    ordered$p_value[[1L]],
    "fit manifest P value after JSON round trip"
  )

  window_manifest <- new_window_run_manifest(
    "analysis-1",
    "chr1_0_2000000",
    ordered,
    c(dosage = "dosage-hash"),
    checkpointed_susie_settings()
  ) |>
    advance_window_run_manifest(committed)
  window_manifest_path <- window_checkpoint_paths(
    "analysis-1",
    "chr1_0_2000000"
  )$window_manifest
  upload_window_run_manifest(store, window_manifest_path, window_manifest)
  hydrated_window_manifest <- read_window_run_manifest_if_present(
    store,
    window_manifest_path
  )
  expect_true_value(
    is_usable_window_manifest(
      hydrated_window_manifest,
      ordered,
      list(
        analysis_id = "analysis-1",
        window_id = "chr1_0_2000000"
      )
    ),
    "window manifest after a very small P value JSON round trip"
  )
})

run_named_test("phenotype commit rejects a checksum mismatch", {
  record <- make_ordered_manifest(1L)[1L, ]
  commit_inputs <- make_commit_inputs(record)
  changed_fit <- fake_fit
  changed_fit$mutation_marker <- TRUE
  saveRDS(changed_fit, commit_inputs$local_artifacts$fit_rds)
  calls <- new.env(parent = emptyenv())
  calls$count <- 0L
  mock_store <- list(
    upload = function(local_path, relative_path) {
      calls$count <- calls$count + 1L
      invisible(TRUE)
    },
    object_uri = function(relative_path) paste0("gs://test-root/", relative_path)
  )

  expect_error_message(
    commit_phenotype_checkpoint(
      mock_store,
      commit_inputs$paths,
      commit_inputs$local_artifacts,
      commit_inputs$fit_manifest
    ),
    "Checksum mismatch"
  )
  expect_identical_value(calls$count, 0L, "uploads after checksum mismatch")
})

run_named_test("phenotype commit validates every Parquet schema before upload", {
  record <- make_ordered_manifest(1L)[1L, ]
  commit_inputs <- make_commit_inputs(record)
  arrow::write_parquet(
    tibble::tibble(wrong_schema = "invalid"),
    commit_inputs$local_artifacts$credible_sets
  )
  commit_inputs$fit_manifest$payloads$credible_sets <- payload_record(
    commit_inputs$local_artifacts$credible_sets,
    commit_inputs$paths$credible_sets
  )
  calls <- new.env(parent = emptyenv())
  calls$count <- 0L
  mock_store <- list(
    upload = function(local_path, relative_path) {
      calls$count <- calls$count + 1L
      invisible(TRUE)
    },
    object_uri = function(relative_path) paste0("gs://test-root/", relative_path)
  )

  expect_error_message(
    commit_phenotype_checkpoint(
      mock_store,
      commit_inputs$paths,
      commit_inputs$local_artifacts,
      commit_inputs$fit_manifest
    ),
    "table columns do not match"
  )
  expect_identical_value(calls$count, 0L, "uploads after Parquet schema failure")
})

run_named_test("phenotype commit enforces every canonical fixed path", {
  record <- make_ordered_manifest(1L)[1L, ]
  path_fields <- c("fit_manifest", "credible_sets", "lbf_variable", "full_susie")

  purrr::walk(path_fields, function(path_field) {
    commit_inputs <- make_commit_inputs(record)
    wrong_path <- file.path("wrong-prefix", basename(commit_inputs$paths[[path_field]]))
    commit_inputs$paths[[path_field]] <- wrong_path
    if (!identical(path_field, "fit_manifest")) {
      commit_inputs$fit_manifest$payloads[[path_field]]$path <- wrong_path
    }
    mock_store <- list(
      upload = function(local_path, relative_path) invisible(TRUE),
      object_uri = function(relative_path) paste0("gs://test-root/", relative_path)
    )

    expect_error_message(
      commit_phenotype_checkpoint(
        mock_store,
        commit_inputs$paths,
        commit_inputs$local_artifacts,
        commit_inputs$fit_manifest
      ),
      "canonical checkpoint path"
    )
  })
})

run_named_test("fit manifest requires the full scientific provenance schema", {
  record <- make_ordered_manifest(1L)[1L, ]
  required_fields <- c(
    "p_value",
    "input_hashes",
    "settings",
    "container_digest",
    "package_versions",
    "source_hashes",
    "runtime",
    "counts",
    "qc",
    "started_at",
    "completed_at"
  )

  purrr::walk(required_fields, function(required_field) {
    commit_inputs <- make_commit_inputs(record)
    commit_inputs$fit_manifest[[required_field]] <- NULL
    mock_store <- list(
      upload = function(local_path, relative_path) invisible(TRUE),
      object_uri = function(relative_path) paste0("gs://test-root/", relative_path)
    )
    expect_error_message(
      commit_phenotype_checkpoint(
        mock_store,
        commit_inputs$paths,
        commit_inputs$local_artifacts,
        commit_inputs$fit_manifest
      ),
      paste0("missing required fields: ", required_field)
    )
  })
})

run_named_test("window manifest helpers advance and record a failure", {
  ordered <- make_ordered_manifest(2L)
  window_manifest <- new_window_run_manifest(
    "analysis-1",
    "chr1_0_2000000",
    ordered,
    c(dosage = "dosage-hash"),
    checkpointed_susie_settings()
  )
  expect_identical_value(window_manifest$last_committed_index, -1L, "initial cursor")
  expect_identical_value(window_manifest$status, "RUNNING", "initial status")

  fit_manifest <- list(
    analysis_id = "analysis-1",
    window_id = "chr1_0_2000000",
    phenotype_id = ordered$phenotype_id[[1L]],
    phenotype_key = ordered$phenotype_key[[1L]],
    modality = ordered$modality[[1L]],
    processing_index = 0L,
    p_value = ordered$p_value[[1L]],
    status = "CONVERGED",
    fit_manifest_path = phenotype_checkpoint_paths(
      "analysis-1",
      "chr1_0_2000000",
      ordered$phenotype_key[[1L]],
      "hash"
    )$fit_manifest
  )
  advanced <- advance_window_run_manifest(window_manifest, fit_manifest)
  expect_identical_value(advanced$last_committed_index, 0L, "advanced cursor")
  expect_identical_value(length(advanced$committed), 1L, "committed record count")

  failed <- fail_window_run_manifest(advanced, simpleError("synthetic failure"))
  expect_identical_value(failed$status, "FAILED", "failed status")
  expect_identical_value(failed$failure$message, "synthetic failure", "failure message")
  expect_identical_value(failed$failure$processing_index, 1L, "failure processing index")
  expect_identical_value(
    failed$failure$phenotype_id,
    ordered$phenotype_id[[2L]],
    "failure phenotype ID"
  )
})

run_named_test("payload recovery keeps a valid recovered-prefix inventory", {
  hydrated <- list(
    recovered_prefix_last_committed_index = 1L,
    last_committed_index = 2L,
    status = "COMPLETE",
    committed = purrr::map(0:2, ~ list(processing_index = .x))
  )
  trimmed_hydrated <- trim_window_manifest_boundary(hydrated, 2L)
  expect_identical_value(
    trimmed_hydrated$recovered_prefix_last_committed_index,
    -1L,
    "hydrated recovery prefix"
  )
  expect_identical_value(
    purrr::map_int(trimmed_hydrated$committed, "processing_index"),
    0:1,
    "hydrated recovery inventory"
  )

  implicit_prefix <- list(
    recovered_prefix_last_committed_index = 1L,
    last_committed_index = 2L,
    status = "RUNNING",
    committed = list(list(processing_index = 2L))
  )
  trimmed_prefix <- trim_window_manifest_boundary(implicit_prefix, 1L)
  expect_identical_value(
    trimmed_prefix$recovered_prefix_last_committed_index,
    0L,
    "implicit recovery prefix"
  )
  expect_identical_value(
    length(trimmed_prefix$committed),
    0L,
    "implicit recovery inventory"
  )
})

run_named_test("window completion uses COMPLETE and invalid statuses are rejected", {
  ordered <- make_ordered_manifest(1L)
  window_manifest <- new_window_run_manifest(
    "analysis-1",
    "chr1_0_2000000",
    ordered,
    c(dosage = paste(rep("a", 64L), collapse = "")),
    checkpointed_susie_settings()
  )
  window_manifest <- advance_window_run_manifest(
    window_manifest,
    list(
      analysis_id = "analysis-1",
      window_id = "chr1_0_2000000",
      phenotype_id = ordered$phenotype_id[[1L]],
      phenotype_key = ordered$phenotype_key[[1L]],
      modality = ordered$modality[[1L]],
      processing_index = 0L,
      p_value = ordered$p_value[[1L]],
      status = "CONVERGED"
    )
  )
  final_record <- list(
    path = "analysis-1/chr1_0_2000000/output",
    sha256 = paste(rep("b", 64L), collapse = ""),
    bytes = 1,
    uri = "gs://test-root/analysis-1/chr1_0_2000000/output"
  )
  final_outputs <- stats::setNames(
    rep(list(final_record), 4L),
    c("fit_index", "credible_sets", "lbf_variable", "full_susie")
  )
  completed <- complete_window_run_manifest(window_manifest, final_outputs)
  expect_identical_value(completed$status, "COMPLETE", "terminal window status")

  purrr::walk(c("COMPLETED", "PENDING", ""), function(invalid_status) {
    invalid_cursor <- new_window_run_manifest(
      "analysis-1",
      "chr1_0_2000000",
      ordered,
      c(dosage = paste(rep("a", 64L), collapse = "")),
      checkpointed_susie_settings()
    )
    invalid_cursor$status <- invalid_status
    resume <- resolve_resume_boundary(
      new_checkpoint_store(tempfile("invalid-window-status-")),
      ordered,
      invalid_cursor
    )
    expect_identical_value(
      resume$window_manifest,
      NULL,
      paste("invalid status", invalid_status)
    )
  })
})

run_named_test("failure manifest uses null identity after the phenotype inventory", {
  ordered <- make_ordered_manifest(1L)
  window_manifest <- new_window_run_manifest(
    "analysis-1",
    "chr1_0_2000000",
    ordered,
    c(dosage = "dosage-hash"),
    checkpointed_susie_settings()
  )
  window_manifest$last_committed_index <- 0L
  failed <- fail_window_run_manifest(window_manifest, simpleError("finalization failure"))
  expect_identical_value(failed$failure$processing_index, NULL, "out-of-range failure index")
  expect_identical_value(failed$failure$phenotype_id, NULL, "out-of-range phenotype ID")
})

run_named_test("corrupt window cursor JSON enters fallback recovery", {
  ordered <- make_ordered_manifest(1L)
  store <- new_checkpoint_store(tempfile("checkpoint-cursor-json-"))
  base_manifest <- new_window_run_manifest(
    "analysis-1",
    "chr1_0_2000000",
    ordered,
    c(dosage = "dosage-hash"),
    checkpointed_susie_settings()
  )

  malformed_key <- base_manifest
  malformed_key$phenotypes[[1L]]$phenotype_key <- c("key_0001", "extra")
  extra_commit <- base_manifest
  extra_commit$committed <- list(list(processing_index = 0L))
  fractional_schema <- base_manifest
  fractional_schema$schema_version <- 1.5
  fractional_index <- base_manifest
  fractional_index$last_committed_index <- -0.5
  fractional_index$committed <- list(list(processing_index = 0L))

  corrupt_manifests <- list(
    malformed_key = malformed_key,
    extra_commit = extra_commit,
    fractional_schema = fractional_schema,
    fractional_index = fractional_index
  )
  purrr::iwalk(corrupt_manifests, function(corrupt_manifest, label) {
    resume <- resolve_resume_boundary(store, ordered, corrupt_manifest)
    expect_identical_value(resume$window_manifest, NULL, paste(label, "fallback cursor"))
    expect_identical_value(resume$next_index, 0L, paste(label, "fallback index"))
  })
})

run_named_test("fractional committed index invalidates the cursor boundary", {
  ordered <- make_ordered_manifest(1L)
  store <- new_checkpoint_store(tempfile("checkpoint-fractional-commit-"))
  commit_inputs <- make_commit_inputs(ordered[1L, ])
  committed <- commit_phenotype_checkpoint(
    store,
    commit_inputs$paths,
    commit_inputs$local_artifacts,
    commit_inputs$fit_manifest
  )
  window_manifest <- new_window_run_manifest(
    "analysis-1",
    "chr1_0_2000000",
    ordered,
    c(dosage = "dosage-hash"),
    checkpointed_susie_settings()
  )
  window_manifest <- advance_window_run_manifest(window_manifest, committed)
  window_manifest$committed[[1L]]$processing_index <- 0.5

  resume <- resolve_resume_boundary(store, ordered, window_manifest)
  expect_identical_value(resume$last_committed_index, 0L, "recovered fixed boundary")
  expect_identical_value(resume$next_index, 1L, "recovered fixed next index")
  expect_identical_value(resume$window_manifest, NULL, "fractional cursor fallback")
})

run_named_test("ordered manifest rejects coordinated cursor identity corruption", {
  ordered <- make_ordered_manifest(1L)
  store <- new_checkpoint_store(tempfile("checkpoint-coordinated-corruption-"))
  commit_inputs <- make_commit_inputs(ordered[1L, ])
  committed <- commit_phenotype_checkpoint(
    store,
    commit_inputs$paths,
    commit_inputs$local_artifacts,
    commit_inputs$fit_manifest
  )
  window_manifest <- new_window_run_manifest(
    "analysis-1",
    "chr1_0_2000000",
    ordered,
    c(dosage = "dosage-hash"),
    checkpointed_susie_settings()
  )
  window_manifest <- advance_window_run_manifest(window_manifest, committed)
  identity_corruption <- window_manifest
  identity_corruption$phenotypes[[1L]]$phenotype_id <- "coordinated-corrupt-id"
  identity_corruption$phenotypes[[1L]]$modality <- "coordinated-corrupt-modality"
  identity_corruption$committed[[1L]]$phenotype_id <- "coordinated-corrupt-id"
  identity_corruption$committed[[1L]]$modality <- "coordinated-corrupt-modality"
  p_value_corruption <- window_manifest
  p_value_corruption$phenotypes[[1L]]$p_value <- 0.75
  p_value_corruption$committed[[1L]]$p_value <- 0.75

  corrupt_cursors <- list(
    coordinated_identity = identity_corruption,
    coordinated_p_value = p_value_corruption
  )
  purrr::iwalk(corrupt_cursors, function(corrupt_cursor, label) {
    resume <- resolve_resume_boundary(store, ordered, corrupt_cursor)
    expect_identical_value(
      resume$last_committed_index,
      0L,
      paste(label, "recovered boundary")
    )
    expect_identical_value(resume$next_index, 1L, paste(label, "recovered next index"))
    expect_identical_value(resume$window_manifest, NULL, paste(label, "fallback cursor"))
  })
})

run_named_test("GCS download failure propagates from resume", {
  ordered <- make_ordered_manifest(1L)
  fake_gsutil <- tempfile("fake-gsutil-")
  writeLines(
    c(
      "#!/bin/sh",
      "if [ \"$2\" = \"stat\" ]; then",
      "  exit 0",
      "fi",
      "exit 17"
    ),
    fake_gsutil
  )
  Sys.chmod(fake_gsutil, mode = "0755")
  store <- new_checkpoint_store("gs://test-bucket/checkpoints", gsutil = fake_gsutil)
  expected_uri <- paste0(
    "gs://test-bucket/checkpoints/",
    "analysis-1/chr1_0_2000000/phenotypes/key_0001/fit_manifest.json"
  )

  expect_error_message(
    resolve_resume_boundary(store, ordered),
    paste0("GCS download failed for: ", expected_uri)
  )
})

run_named_test("fallback recovery uses logarithmic exact existence checks", {
  ordered <- make_ordered_manifest(1024L)
  store_root <- tempfile("checkpoint-prefix-")
  store <- new_checkpoint_store(store_root)

  purrr::walk(0:730, function(index) {
    record <- ordered[index + 1L, ]
    manifest_path <- phenotype_checkpoint_paths(
      record$analysis_id[[1L]],
      record$window_id[[1L]],
      record$phenotype_key[[1L]]
    )$fit_manifest
    destination <- file.path(store_root, manifest_path)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    writeLines("{}", destination)
  })

  boundary_record <- ordered[731L, ]
  skipped_inputs <- make_skipped_commit_inputs(boundary_record)
  commit_phenotype_checkpoint(
    store,
    skipped_inputs$paths,
    skipped_inputs$local_artifacts,
    skipped_inputs$fit_manifest
  )

  check_counter <- new.env(parent = emptyenv())
  check_counter$count <- 0L
  counting_store <- store
  counting_store$object_exists <- function(relative_path) {
    check_counter$count <- check_counter$count + 1L
    store$object_exists(relative_path)
  }

  resume <- resolve_resume_boundary(counting_store, ordered)
  expect_identical_value(resume$last_committed_index, 730L, "recovered boundary")
  expect_identical_value(resume$next_index, 731L, "recovered next index")
  expect_identical_value(resume$recovered, TRUE, "fallback recovery flag")
  expect_identical_value(resume$window_manifest, NULL, "fallback window manifest")
  expect_true_value(check_counter$count <= 12L, "logarithmic existence-check count")
  cat("INFO boundary existence checks: ", check_counter$count, "\n", sep = "")
})

run_named_test("corrupt boundary is returned for recomputation", {
  ordered <- make_ordered_manifest(1L)
  store_root <- tempfile("checkpoint-corrupt-")
  store <- new_checkpoint_store(store_root)
  commit_inputs <- make_commit_inputs(ordered[1L, ])
  committed <- commit_phenotype_checkpoint(
    store,
    commit_inputs$paths,
    commit_inputs$local_artifacts,
    commit_inputs$fit_manifest
  )
  committed$payloads$susie_fit$sha256 <- paste(rep("f", 64L), collapse = "")
  jsonlite::write_json(
    committed,
    file.path(store_root, commit_inputs$paths$fit_manifest),
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )

  window_manifest <- new_window_run_manifest(
    "analysis-1",
    "chr1_0_2000000",
    ordered,
    c(dosage = "dosage-hash"),
    checkpointed_susie_settings()
  )
  window_manifest <- advance_window_run_manifest(window_manifest, committed)

  resume <- resolve_resume_boundary(store, ordered, window_manifest)
  expect_identical_value(resume$last_committed_index, -1L, "corrupt prior boundary")
  expect_identical_value(resume$next_index, 0L, "corrupt recomputation index")
  expect_identical_value(
    resume$window_manifest$last_committed_index,
    -1L,
    "repaired corrupt cursor"
  )
  expect_identical_value(
    length(resume$window_manifest$recovery_history),
    1L,
    "corruption recovery history length"
  )
  recovery <- resume$window_manifest$recovery_history[[1L]]
  expect_identical_value(recovery$processing_index, 0L, "recovery index")
  expect_identical_value(recovery$phenotype_id, ordered$phenotype_id[[1L]], "recovery phenotype")
  expect_identical_value(recovery$phenotype_key, ordered$phenotype_key[[1L]], "recovery key")
  expect_identical_value(
    recovery$fit_manifest_path,
    commit_inputs$paths$fit_manifest,
    "recovery manifest path"
  )
  expect_true_value(nzchar(recovery$reason), "recovery reason")
})

run_named_test("fallback rejects a corrupt nonconverged RDS", {
  ordered <- make_ordered_manifest(1L)
  store_root <- tempfile("checkpoint-nonconverged-")
  store <- new_checkpoint_store(store_root)
  nonconverged_fit <- fake_fit
  nonconverged_fit$converged <- FALSE
  commit_inputs <- make_commit_inputs(ordered[1L, ], nonconverged_fit)
  commit_phenotype_checkpoint(
    store,
    commit_inputs$paths,
    commit_inputs$local_artifacts,
    commit_inputs$fit_manifest
  )
  writeBin(charToRaw("corrupt RDS"), file.path(store_root, commit_inputs$paths$fit_rds))

  resume <- resolve_resume_boundary(store, ordered)
  expect_identical_value(resume$last_committed_index, -1L, "nonconverged corrupt boundary")
  expect_identical_value(resume$next_index, 0L, "nonconverged recomputation index")
})

run_named_test("fallback rejects an unexpected skipped reason", {
  ordered <- make_ordered_manifest(1L)
  store_root <- tempfile("checkpoint-skipped-reason-")
  store <- new_checkpoint_store(store_root)
  record <- ordered[1L, ]
  skipped_inputs <- make_skipped_commit_inputs(record)
  committed <- commit_phenotype_checkpoint(
    store,
    skipped_inputs$paths,
    skipped_inputs$local_artifacts,
    skipped_inputs$fit_manifest
  )
  committed$exclusion_reason <- "UNEXPECTED_REASON"
  jsonlite::write_json(
    committed,
    file.path(store_root, skipped_inputs$paths$fit_manifest),
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )

  resume <- resolve_resume_boundary(store, ordered)
  expect_identical_value(resume$last_committed_index, -1L, "invalid skipped boundary")
  expect_identical_value(resume$next_index, 0L, "invalid skipped recomputation index")
})

run_named_test("resume adopts one valid commit after the saved cursor", {
  ordered <- make_ordered_manifest(12L)
  store_root <- tempfile("checkpoint-cursor-")
  store <- new_checkpoint_store(store_root)

  commit_9 <- make_commit_inputs(ordered[10L, ])
  committed_9 <- commit_phenotype_checkpoint(
    store,
    commit_9$paths,
    commit_9$local_artifacts,
    commit_9$fit_manifest
  )
  commit_10 <- make_commit_inputs(ordered[11L, ])
  commit_phenotype_checkpoint(
    store,
    commit_10$paths,
    commit_10$local_artifacts,
    commit_10$fit_manifest
  )

  window_manifest <- new_window_run_manifest(
    "analysis-1",
    "chr1_0_2000000",
    ordered,
    c(dosage = "dosage-hash"),
    checkpointed_susie_settings()
  )
  purrr::walk(0:8, function(index) {
    record <- ordered[index + 1L, ]
    window_manifest <<- advance_window_run_manifest(
      window_manifest,
      list(
        analysis_id = record$analysis_id[[1L]],
        window_id = record$window_id[[1L]],
        phenotype_id = record$phenotype_id[[1L]],
        phenotype_key = record$phenotype_key[[1L]],
        modality = record$modality[[1L]],
        processing_index = record$processing_index[[1L]],
        p_value = record$p_value[[1L]],
        status = "CONVERGED",
        fit_manifest_path = phenotype_checkpoint_paths(
          record$analysis_id[[1L]],
          record$window_id[[1L]],
          record$phenotype_key[[1L]],
          "prior"
        )$fit_manifest
      )
    )
  })
  window_manifest <- advance_window_run_manifest(window_manifest, committed_9)

  resume <- resolve_resume_boundary(store, ordered, window_manifest)
  expect_identical_value(resume$last_committed_index, 10L, "adopted boundary")
  expect_identical_value(resume$next_index, 11L, "adopted next index")
  expect_identical_value(resume$recovered, TRUE, "adopted commit flag")
  expect_identical_value(resume$window_manifest$last_committed_index, 10L, "updated cursor")
})

run_named_test("resume records an invalid exact next manifest before recomputation", {
  ordered <- make_ordered_manifest(2L)
  store_root <- tempfile("checkpoint-invalid-next-")
  on.exit(unlink(store_root, recursive = TRUE), add = TRUE)
  store <- new_checkpoint_store(store_root)

  first_inputs <- make_commit_inputs(ordered[1L, ])
  first_committed <- commit_phenotype_checkpoint(
    store,
    first_inputs$paths,
    first_inputs$local_artifacts,
    first_inputs$fit_manifest
  )
  second_inputs <- make_commit_inputs(ordered[2L, ])
  second_committed <- commit_phenotype_checkpoint(
    store,
    second_inputs$paths,
    second_inputs$local_artifacts,
    second_inputs$fit_manifest
  )
  second_committed$payloads$susie_fit$sha256 <- paste(rep("f", 64L), collapse = "")
  jsonlite::write_json(
    second_committed,
    file.path(store_root, second_inputs$paths$fit_manifest),
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null",
    digits = NA
  )

  window_manifest <- new_window_run_manifest(
    "analysis-1",
    "chr1_0_2000000",
    ordered,
    c(dosage = "dosage-hash"),
    checkpointed_susie_settings()
  ) |>
    advance_window_run_manifest(first_committed)
  resume <- resolve_resume_boundary(store, ordered, window_manifest)

  expect_identical_value(resume$last_committed_index, 0L, "retained valid boundary")
  expect_identical_value(resume$next_index, 1L, "invalid next recomputation index")
  expect_identical_value(resume$recovered, TRUE, "audited recovery state")
  expect_identical_value(
    length(resume$window_manifest$recovery_history),
    1L,
    "invalid next recovery count"
  )
  recovery <- resume$window_manifest$recovery_history[[1L]]
  expect_identical_value(recovery$processing_index, 1L, "invalid next recovery index")
  expect_identical_value(
    recovery$fit_manifest_path,
    second_inputs$paths$fit_manifest,
    "invalid next recovery path"
  )
  expect_true_value(nzchar(recovery$reason), "invalid next recovery reason")
  expect_true_value(
    store$object_exists(second_inputs$paths$fit_manifest),
    "invalid next manifest remains durable"
  )
})

run_named_test("resume compares name-stable current provenance with cursor and boundary", {
  expected_input_hashes <- c(
    dosage = paste(rep("a", 64L), collapse = ""),
    phenotypes = paste(rep("b", 64L), collapse = "")
  )
  expected_settings <- checkpointed_susie_settings()

  make_provenance_checkpoint <- function(label) {
    ordered <- make_ordered_manifest(1L)
    store_root <- tempfile(paste0("checkpoint-provenance-", label, "-"))
    store <- new_checkpoint_store(store_root)
    commit_inputs <- make_commit_inputs(ordered[1L, ])
    committed <- commit_phenotype_checkpoint(
      store,
      commit_inputs$paths,
      commit_inputs$local_artifacts,
      commit_inputs$fit_manifest
    )
    cursor <- new_window_run_manifest(
      "analysis-1",
      "chr1_0_2000000",
      ordered,
      expected_input_hashes,
      expected_settings
    ) |>
      advance_window_run_manifest(committed)
    list(
      ordered = ordered,
      store = store,
      store_root = store_root,
      committed = committed,
      cursor = cursor,
      fit_manifest_path = commit_inputs$paths$fit_manifest
    )
  }

  unchanged <- make_provenance_checkpoint("unchanged")
  on.exit(unlink(unchanged$store_root, recursive = TRUE), add = TRUE)
  stable <- resolve_resume_boundary(
    unchanged$store,
    unchanged$ordered,
    unchanged$cursor,
    expected_input_hashes = rev(expected_input_hashes),
    expected_settings = expected_settings[rev(names(expected_settings))]
  )
  expect_identical_value(stable$recovered, FALSE, "name-stable provenance")
  expect_identical_value(stable$next_index, 1L, "name-stable next index")

  purrr::walk(c("input_hashes", "settings"), function(field) {
    window_only <- make_provenance_checkpoint(paste0("window-", field))
    on.exit(unlink(window_only$store_root, recursive = TRUE), add = TRUE)
    if (identical(field, "input_hashes")) {
      window_only$cursor$input_hashes$dosage <- paste(rep("c", 64L), collapse = "")
    } else {
      window_only$cursor$settings$L <- 11L
    }
    recovered <- resolve_resume_boundary(
      window_only$store,
      window_only$ordered,
      window_only$cursor,
      expected_input_hashes = expected_input_hashes,
      expected_settings = expected_settings
    )
    expect_identical_value(recovered$window_manifest, NULL, paste(field, "cursor fallback"))
    expect_identical_value(recovered$next_index, 1L, paste(field, "fixed checkpoint reuse"))

    coordinated <- make_provenance_checkpoint(paste0("coordinated-", field))
    on.exit(unlink(coordinated$store_root, recursive = TRUE), add = TRUE)
    if (identical(field, "input_hashes")) {
      coordinated$cursor$input_hashes$dosage <- paste(rep("c", 64L), collapse = "")
      coordinated$committed$input_hashes$dosage <- paste(rep("c", 64L), collapse = "")
    } else {
      coordinated$cursor$settings$L <- 11L
      coordinated$committed$settings$L <- 11L
    }
    jsonlite::write_json(
      coordinated$committed,
      file.path(coordinated$store_root, coordinated$fit_manifest_path),
      auto_unbox = TRUE,
      pretty = TRUE,
      null = "null",
      digits = NA
    )
    recomputed <- resolve_resume_boundary(
      coordinated$store,
      coordinated$ordered,
      coordinated$cursor,
      expected_input_hashes = expected_input_hashes,
      expected_settings = expected_settings
    )
    expect_identical_value(recomputed$last_committed_index, -1L, paste(field, "recompute boundary"))
    expect_identical_value(recomputed$next_index, 0L, paste(field, "recompute index"))
  })
})
