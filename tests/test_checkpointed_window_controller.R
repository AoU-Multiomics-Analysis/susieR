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
    wrapper_path = controller_script
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
    sha256_file(config$wrapper_path)
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
    checkpointed_susie_settings()
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
      "checkpoint_root", "output_dir", "wrapper_path"
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
