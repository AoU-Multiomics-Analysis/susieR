source("tests/test_helpers.R")
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
  purrr::iwalk(
    paths[c("credible_sets", "lbf_variable", "full_susie")],
    ~ arrow::write_parquet(tibble::tibble(artifact = .y), .x)
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
    status = if (isTRUE(fit$converged)) "CONVERGED" else "NONCONVERGED",
    converged = fit$converged,
    exclusion_reason = NULL,
    payloads = list(
      susie_fit = payload_record(local_artifacts$fit_rds, paths$fit_rds),
      credible_sets = payload_record(local_artifacts$credible_sets, paths$credible_sets),
      lbf_variable = payload_record(local_artifacts$lbf_variable, paths$lbf_variable),
      full_susie = payload_record(local_artifacts$full_susie, paths$full_susie)
    )
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

run_named_test("phenotype commit rejects a checksum mismatch", {
  record <- make_ordered_manifest(1L)[1L, ]
  commit_inputs <- make_commit_inputs(record)
  commit_inputs$fit_manifest$payloads$susie_fit$sha256 <- paste(rep("0", 64L), collapse = "")
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
  boundary_path <- phenotype_checkpoint_paths(
    boundary_record$analysis_id[[1L]],
    boundary_record$window_id[[1L]],
    boundary_record$phenotype_key[[1L]]
  )$fit_manifest
  jsonlite::write_json(
    list(
      schema_version = 1L,
      analysis_id = boundary_record$analysis_id[[1L]],
      window_id = boundary_record$window_id[[1L]],
      phenotype_id = boundary_record$phenotype_id[[1L]],
      phenotype_key = boundary_record$phenotype_key[[1L]],
      modality = boundary_record$modality[[1L]],
      processing_index = boundary_record$processing_index[[1L]],
      status = "SKIPPED",
      exclusion_reason = "NO_USABLE_VARIANTS",
      converged = NULL,
      payloads = list()
    ),
    file.path(store_root, boundary_path),
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
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
  manifest_path <- phenotype_checkpoint_paths(
    record$analysis_id[[1L]],
    record$window_id[[1L]],
    record$phenotype_key[[1L]]
  )$fit_manifest
  destination <- file.path(store_root, manifest_path)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(
    list(
      schema_version = 1L,
      analysis_id = record$analysis_id[[1L]],
      window_id = record$window_id[[1L]],
      phenotype_id = record$phenotype_id[[1L]],
      phenotype_key = record$phenotype_key[[1L]],
      modality = record$modality[[1L]],
      processing_index = record$processing_index[[1L]],
      status = "SKIPPED",
      exclusion_reason = "UNEXPECTED_REASON",
      converged = NULL,
      payloads = list()
    ),
    destination,
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
