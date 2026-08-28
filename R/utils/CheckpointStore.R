# Exact-path checkpoint storage and resume helpers.

checkpoint_schema_version <- 1L

checkpoint_terminal_statuses <- c("CONVERGED", "NONCONVERGED", "SKIPPED")

checkpoint_window_statuses <- c("RUNNING", "FAILED", "COMPLETE")

checkpoint_skip_reasons <- c(
  "NO_USABLE_VARIANTS",
  "ZERO_PHENOTYPE_VARIANCE",
  "NO_ALTERNATIVE_ALLELE",
  "TOO_FEW_ALIGNED_SAMPLES"
)

new_checkpoint_store_error <- function(message) {
  structure(
    list(message = message, call = NULL),
    class = c("checkpoint_store_error", "error", "condition")
  )
}

stop_checkpoint_store_error <- function(message) {
  stop(new_checkpoint_store_error(message))
}

new_checkpoint_store <- function(root, gsutil = "gsutil") {
  is_gcs <- startsWith(root, "gs://")
  clean_root <- sub("/+$", "", root)
  object_uri <- function(relative_path) paste0(clean_root, "/", relative_path)
  capture_gsutil <- function(arguments, operation, uri) {
    output <- tryCatch(
      suppressWarnings(system2(gsutil, arguments, stdout = TRUE, stderr = TRUE)),
      error = function(condition) {
        stop_checkpoint_store_error(
          paste0("GCS ", operation, " failed for: ", uri)
        )
      }
    )
    list(
      status = as.integer(attr(output, "status") %||% 0L),
      output = paste(output, collapse = "\n")
    )
  }
  run_gsutil <- function(arguments, operation, uri) {
    result <- capture_gsutil(arguments, operation, uri)
    if (!identical(result$status, 0L)) {
      stop_checkpoint_store_error(paste0("GCS ", operation, " failed for: ", uri))
    }
    invisible(TRUE)
  }

  list(
    root = clean_root,
    is_gcs = is_gcs,
    object_exists = function(relative_path) {
      if (is_gcs) {
        uri <- object_uri(relative_path)
        result <- capture_gsutil(c("-q", "stat", uri), "stat", uri)
        if (identical(result$status, 0L)) {
          return(TRUE)
        }
        not_found_pattern <- paste(
          c(
            "No URLs matched",
            "matched no objects",
            "No such object",
            "NotFoundException",
            "(^|[^0-9])404([^0-9]|$)"
          ),
          collapse = "|"
        )
        if (grepl(not_found_pattern, result$output, ignore.case = TRUE)) {
          return(FALSE)
        }
        stop_checkpoint_store_error(paste0("GCS stat failed for: ", uri))
      }
      file.exists(file.path(clean_root, relative_path))
    },
    upload = function(local_path, relative_path) {
      if (is_gcs) {
        return(run_gsutil(c("-q", "cp", local_path, object_uri(relative_path)), "upload", object_uri(relative_path)))
      }
      destination <- file.path(clean_root, relative_path)
      dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
      if (!file.copy(local_path, destination, overwrite = TRUE)) {
        stop_checkpoint_store_error(paste0("Filesystem upload failed for: ", destination))
      }
      invisible(TRUE)
    },
    download = function(relative_path, local_path) {
      tryCatch(
        dir.create(dirname(local_path), recursive = TRUE, showWarnings = FALSE),
        error = function(condition) {
          stop_checkpoint_store_error(
            paste0("Checkpoint download setup failed for: ", local_path)
          )
        }
      )
      if (is_gcs) {
        return(run_gsutil(c("-q", "cp", object_uri(relative_path), local_path), "download", object_uri(relative_path)))
      }
      source_path <- file.path(clean_root, relative_path)
      copy_succeeded <- tryCatch(
        file.copy(source_path, local_path, overwrite = TRUE),
        error = function(condition) FALSE
      )
      if (!copy_succeeded) {
        stop_checkpoint_store_error(paste0("Filesystem download failed for: ", source_path))
      }
      invisible(TRUE)
    },
    object_uri = object_uri
  )
}

read_window_run_manifest_if_present <- function(store, relative_path) {
  if (!store$object_exists(relative_path)) {
    return(NULL)
  }
  local_manifest <- tempfile("window-run-manifest-", fileext = ".json")
  on.exit(unlink(local_manifest), add = TRUE)
  store$download(relative_path, local_manifest)
  window_manifest <- tryCatch(
    jsonlite::read_json(local_manifest, simplifyVector = FALSE),
    error = function(condition) NULL
  )
  if (!is.list(window_manifest)) {
    return(NULL)
  }
  window_manifest
}

upload_window_run_manifest <- function(store, relative_path, window_manifest) {
  if (!is.list(window_manifest)) {
    stop("Window run manifest must be a JSON object.", call. = FALSE)
  }
  local_manifest <- tempfile("window-run-manifest-", fileext = ".json")
  on.exit(unlink(local_manifest), add = TRUE)
  jsonlite::write_json(
    window_manifest,
    local_manifest,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null",
    digits = NA
  )
  store$upload(local_manifest, relative_path)
  TRUE
}

window_checkpoint_paths <- function(analysis_id, window_id) {
  list(
    window_manifest = file.path(analysis_id, window_id, "window_manifest.json"),
    fit_index = file.path(analysis_id, window_id, "window_fit_index.tsv"),
    credible_sets = file.path(analysis_id, window_id, "results", "credible_sets.parquet"),
    lbf_variable = file.path(analysis_id, window_id, "results", "lbf_variable.parquet"),
    full_susie = file.path(analysis_id, window_id, "results", "full_susie.parquet")
  )
}

phenotype_checkpoint_paths <- function(
    analysis_id,
    window_id,
    phenotype_key,
    fit_sha256 = NULL
) {
  phenotype_root <- file.path(
    analysis_id,
    window_id,
    "phenotypes",
    phenotype_key
  )
  fit_rds <- if (is.null(fit_sha256)) {
    character(0)
  } else {
    file.path(phenotype_root, paste0("susie_fit.", fit_sha256, ".rds"))
  }

  list(
    fit_manifest = file.path(phenotype_root, "fit_manifest.json"),
    fit_rds = fit_rds,
    credible_sets = file.path(phenotype_root, "credible_sets.parquet"),
    lbf_variable = file.path(phenotype_root, "lbf_variable.parquet"),
    full_susie = file.path(phenotype_root, "full_susie.parquet")
  )
}

checkpoint_scalar_character <- function(value, label) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop(label, " must be one nonempty string.", call. = FALSE)
  }
  value
}

checkpoint_scalar_index <- function(value, label) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value != floor(value) ||
      value < -.Machine$integer.max || value > .Machine$integer.max) {
    stop(label, " must be one integer.", call. = FALSE)
  }
  as.integer(value)
}

checkpoint_optional_index <- function(value) {
  tryCatch(
    checkpoint_scalar_index(value, "Checkpoint index"),
    error = function(condition) NULL
  )
}

validate_checkpoint_named_strings <- function(value, label, pattern = NULL) {
  if (!is.list(value) || length(value) == 0L || is.null(names(value)) ||
      anyNA(names(value)) || any(!nzchar(names(value))) || anyDuplicated(names(value))) {
    stop(label, " must be a nonempty named object.", call. = FALSE)
  }
  purrr::iwalk(value, function(item, name) {
    checkpoint_scalar_character(item, paste0(label, " ", name))
    if (!is.null(pattern) && !grepl(pattern, item)) {
      stop(label, " ", name, " is invalid.", call. = FALSE)
    }
  })
  invisible(TRUE)
}

validate_checkpoint_counts <- function(counts) {
  required_counts <- c(
    "input_samples",
    "retained_samples",
    "raw_dosage_samples",
    "overlap_samples",
    "phenotype_retained_samples",
    "input_variants",
    "retained_variants"
  )
  if (!is.list(counts)) {
    stop("Checkpoint counts must be an object.", call. = FALSE)
  }
  missing_counts <- setdiff(required_counts, names(counts))
  if (length(missing_counts) > 0L) {
    stop(
      "Checkpoint counts are missing required fields: ",
      paste(missing_counts, collapse = ", "),
      call. = FALSE
    )
  }
  validated_counts <- purrr::map(
    counts[required_counts],
    ~ checkpoint_scalar_index(.x, "Checkpoint count")
  )
  if (any(unlist(validated_counts, use.names = FALSE) < 0L)) {
    stop("Checkpoint counts cannot be negative.", call. = FALSE)
  }
  if (validated_counts$retained_samples > validated_counts$input_samples ||
      validated_counts$overlap_samples > validated_counts$raw_dosage_samples ||
      validated_counts$phenotype_retained_samples > validated_counts$overlap_samples ||
      validated_counts$retained_variants > validated_counts$input_variants) {
    stop("Checkpoint retained counts cannot exceed input counts.", call. = FALSE)
  }
  invisible(TRUE)
}

validate_checkpoint_qc <- function(qc) {
  required <- c(
    "removed_covariates",
    "removed_variant_ids",
    "design_id",
    "sample_mask_id",
    "cache_key"
  )
  if (!is.list(qc) || any(!required %in% names(qc))) {
    stop("Checkpoint QC is missing required fields.", call. = FALSE)
  }
  purrr::walk(c("removed_covariates", "removed_variant_ids"), function(field) {
    values <- unlist(qc[[field]], use.names = FALSE)
    if (!is.character(values) || anyNA(values) || any(!nzchar(values))) {
      if (length(values) > 0L) {
        stop("Checkpoint QC identifier lists must contain nonempty strings.", call. = FALSE)
      }
    }
  })
  purrr::walk(c("design_id", "sample_mask_id", "cache_key"), function(field) {
    value <- checkpoint_scalar_character(qc[[field]], paste0("Checkpoint QC ", field))
    if (!grepl("^[0-9a-f]{64}$", value)) {
      stop("Checkpoint QC hash is invalid: ", field, call. = FALSE)
    }
  })
  invisible(TRUE)
}

validate_checkpoint_runtime <- function(runtime) {
  if (!is.list(runtime)) {
    stop("Checkpoint runtime identity must be an object.", call. = FALSE)
  }
  checkpoint_scalar_character(runtime$runner_image, "Checkpoint runner image")
  checkpoint_scalar_character(
    runtime$base_image_digest,
    "Checkpoint base image digest"
  )
  invisible(TRUE)
}

validate_checkpoint_payload_record <- function(payload, label) {
  if (!is.list(payload)) {
    stop("Checkpoint payload metadata is missing for ", label, ".", call. = FALSE)
  }
  checkpoint_scalar_character(payload$path, paste0("Checkpoint ", label, " path"))
  checksum <- checkpoint_scalar_character(
    payload$sha256,
    paste0("Checkpoint ", label, " SHA-256")
  )
  if (!grepl("^[0-9a-f]{64}$", checksum)) {
    stop("Checkpoint ", label, " SHA-256 is invalid.", call. = FALSE)
  }
  bytes <- suppressWarnings(as.numeric(payload$bytes))
  if (length(bytes) != 1L || is.na(bytes) || !is.finite(bytes) ||
      bytes < 0 || bytes != floor(bytes)) {
    stop("Checkpoint ", label, " byte size is invalid.", call. = FALSE)
  }
  invisible(TRUE)
}

validate_fit_manifest_structure <- function(fit_manifest, record = NULL) {
  if (!is.list(fit_manifest)) {
    stop("Checkpoint fit manifest must be a JSON object.", call. = FALSE)
  }
  required_fields <- c(
    "schema_version",
    "analysis_id",
    "window_id",
    "phenotype_id",
    "phenotype_key",
    "modality",
    "processing_index",
    "p_value",
    "status",
    "converged",
    "exclusion_reason",
    "input_hashes",
    "settings",
    "container_digest",
    "package_versions",
    "source_hashes",
    "runtime",
    "counts",
    "qc",
    "started_at",
    "completed_at",
    "payloads"
  )
  missing_fields <- setdiff(required_fields, names(fit_manifest))
  if (length(missing_fields) > 0L) {
    stop(
      "Checkpoint fit manifest is missing required fields: ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }
  schema_version <- checkpoint_scalar_index(
    fit_manifest$schema_version,
    "Checkpoint schema version"
  )
  if (!identical(schema_version, checkpoint_schema_version)) {
    stop("Checkpoint fit manifest schema version is not supported.", call. = FALSE)
  }

  identity_fields <- c(
    "analysis_id",
    "window_id",
    "phenotype_id",
    "phenotype_key",
    "modality"
  )
  purrr::walk(
    identity_fields,
    ~ checkpoint_scalar_character(
      fit_manifest[[.x]],
      paste0("Checkpoint ", .x)
    )
  )
  processing_index <- checkpoint_scalar_index(
    fit_manifest$processing_index,
    "Checkpoint processing index"
  )
  p_value <- fit_manifest$p_value
  if (!is.numeric(p_value) || length(p_value) != 1L || is.na(p_value) ||
      !is.finite(p_value) || p_value < 0 || p_value > 1) {
    stop("Checkpoint p_value must be finite and between zero and one.", call. = FALSE)
  }
  validate_checkpoint_named_strings(
    fit_manifest$input_hashes,
    "Checkpoint input hashes",
    "^[0-9a-f]{64}$"
  )
  if (!is.list(fit_manifest$settings) || length(fit_manifest$settings) == 0L ||
      is.null(names(fit_manifest$settings)) || any(!nzchar(names(fit_manifest$settings)))) {
    stop("Checkpoint settings must be a nonempty named object.", call. = FALSE)
  }
  checkpoint_scalar_character(
    fit_manifest$container_digest,
    "Checkpoint container digest"
  )
  validate_checkpoint_named_strings(
    fit_manifest$package_versions,
    "Checkpoint package versions"
  )
  validate_checkpoint_named_strings(
    fit_manifest$source_hashes,
    "Checkpoint runtime source hashes",
    "^[0-9a-f]{64}$"
  )
  validate_checkpoint_runtime(fit_manifest$runtime)
  validate_checkpoint_counts(fit_manifest$counts)
  validate_checkpoint_qc(fit_manifest$qc)
  checkpoint_scalar_character(fit_manifest$started_at, "Checkpoint start timestamp")
  checkpoint_scalar_character(fit_manifest$completed_at, "Checkpoint completion timestamp")
  status <- checkpoint_scalar_character(fit_manifest$status, "Checkpoint status")
  if (!status %in% checkpoint_terminal_statuses) {
    stop("Checkpoint status is not terminal.", call. = FALSE)
  }

  if (!is.null(record)) {
    purrr::walk(identity_fields, function(field) {
      expected <- as.character(record[[field]][[1L]])
      if (!identical(as.character(fit_manifest[[field]]), expected)) {
        stop("Checkpoint ", field, " does not match the ordered manifest.", call. = FALSE)
      }
    })
    expected_index <- checkpoint_scalar_index(
      record$processing_index[[1L]],
      "Ordered manifest processing index"
    )
    if (!identical(processing_index, expected_index)) {
      stop("Checkpoint processing index does not match the ordered manifest.", call. = FALSE)
    }
    expected_p_value <- as.numeric(record$p_value[[1L]])
    if (!checkpoint_p_value_matches(fit_manifest$p_value, expected_p_value)) {
      stop("Checkpoint p_value does not match the ordered manifest.", call. = FALSE)
    }
  }

  if (!is.list(fit_manifest$payloads)) {
    stop("Checkpoint payloads must be an object.", call. = FALSE)
  }
  parquet_payloads <- c("credible_sets", "lbf_variable", "full_susie")
  purrr::walk(
    parquet_payloads,
    ~ validate_checkpoint_payload_record(fit_manifest$payloads[[.x]], .x)
  )

  if (identical(status, "SKIPPED")) {
    reason <- checkpoint_scalar_character(
      fit_manifest$exclusion_reason,
      "Checkpoint exclusion reason"
    )
    if (!reason %in% checkpoint_skip_reasons) {
      stop("Checkpoint exclusion reason is not allowed.", call. = FALSE)
    }
    if (!is.null(fit_manifest$payloads$susie_fit)) {
      stop("A skipped checkpoint cannot contain an RDS payload.", call. = FALSE)
    }
  } else {
    if (!is.logical(fit_manifest$converged) ||
        length(fit_manifest$converged) != 1L || is.na(fit_manifest$converged)) {
      stop("Checkpoint convergence must be one logical value.", call. = FALSE)
    }
    expected_convergence <- identical(status, "CONVERGED")
    if (!identical(fit_manifest$converged, expected_convergence)) {
      stop("Checkpoint status and convergence do not match.", call. = FALSE)
    }
    validate_checkpoint_payload_record(fit_manifest$payloads$susie_fit, "susie_fit")
  }

  fit_sha256 <- if (identical(status, "SKIPPED")) {
    NULL
  } else {
    fit_manifest$payloads$susie_fit$sha256
  }
  canonical_paths <- phenotype_checkpoint_paths(
    fit_manifest$analysis_id,
    fit_manifest$window_id,
    fit_manifest$phenotype_key,
    fit_sha256
  )
  payload_mappings <- checkpoint_artifact_mappings(status)
  purrr::iwalk(payload_mappings, function(payload_name, artifact_name) {
    if (!identical(fit_manifest$payloads[[payload_name]]$path, canonical_paths[[artifact_name]])) {
      stop(
        "Checkpoint payload does not use the canonical checkpoint path: ",
        artifact_name,
        call. = FALSE
      )
    }
  })
  invisible(TRUE)
}

checkpoint_artifact_mappings <- function(status) {
  mapping <- c(
    fit_rds = "susie_fit",
    credible_sets = "credible_sets",
    lbf_variable = "lbf_variable",
    full_susie = "full_susie"
  )
  if (identical(status, "SKIPPED")) {
    mapping[names(mapping) != "fit_rds"]
  } else {
    mapping
  }
}

commit_phenotype_checkpoint <- function(store, paths, local_artifacts, fit_manifest) {
  validate_fit_manifest_structure(fit_manifest)
  status <- fit_manifest$status
  artifact_mappings <- checkpoint_artifact_mappings(status)
  fit_sha256 <- if (identical(status, "SKIPPED")) {
    NULL
  } else {
    fit_manifest$payloads$susie_fit$sha256
  }
  canonical_paths <- phenotype_checkpoint_paths(
    fit_manifest$analysis_id,
    fit_manifest$window_id,
    fit_manifest$phenotype_key,
    fit_sha256
  )
  purrr::walk(names(canonical_paths), function(path_name) {
    if (!identical(paths[[path_name]], canonical_paths[[path_name]])) {
      stop(
        "Path is not the canonical checkpoint path: ",
        path_name,
        call. = FALSE
      )
    }
  })

  purrr::iwalk(artifact_mappings, function(payload_name, artifact_name) {
    local_path <- local_artifacts[[artifact_name]]
    relative_path <- paths[[artifact_name]]
    if (!is.character(local_path) || length(local_path) != 1L || !file.exists(local_path)) {
      stop("Local checkpoint artifact is absent: ", artifact_name, call. = FALSE)
    }
    if (!is.character(relative_path) || length(relative_path) != 1L || !nzchar(relative_path)) {
      stop("Checkpoint path is absent: ", artifact_name, call. = FALSE)
    }
    payload <- fit_manifest$payloads[[payload_name]]
    if (!identical(payload$path, relative_path)) {
      stop("Checkpoint payload path mismatch for: ", artifact_name, call. = FALSE)
    }
    if (!identical(sha256_file(local_path), payload$sha256)) {
      stop("Checksum mismatch for checkpoint artifact: ", artifact_name, call. = FALSE)
    }
    if (!identical(as.numeric(unname(file.info(local_path)$size)), as.numeric(payload$bytes))) {
      stop("Byte-size mismatch for checkpoint artifact: ", artifact_name, call. = FALSE)
    }
  })

  parquet_names <- c("credible_sets", "lbf_variable", "full_susie")
  local_tables <- purrr::map(parquet_names, function(payload_name) {
    artifact_name <- names(artifact_mappings)[
      match(payload_name, unname(artifact_mappings))
    ]
    arrow::read_parquet(
      local_artifacts[[artifact_name]],
      as_data_frame = TRUE
    )
  }) |>
    stats::setNames(parquet_names)
  checkpointed_validate_susie_tables(local_tables)

  if (!identical(status, "SKIPPED")) {
    fit <- tryCatch(
      readRDS(local_artifacts$fit_rds),
      error = function(condition) {
        stop("Checkpoint RDS validation failed: ", conditionMessage(condition), call. = FALSE)
      }
    )
    validate_checkpointed_susie_fit(fit)
    if (!identical(fit$converged, fit_manifest$converged)) {
      stop("Checkpoint fit convergence does not match its manifest.", call. = FALSE)
    }
    expected_fit_path <- phenotype_checkpoint_paths(
      fit_manifest$analysis_id,
      fit_manifest$window_id,
      fit_manifest$phenotype_key,
      fit_manifest$payloads$susie_fit$sha256
    )$fit_rds
    if (!identical(paths$fit_rds, expected_fit_path)) {
      stop("Checkpoint RDS path does not match its checksum.", call. = FALSE)
    }
  }

  purrr::iwalk(artifact_mappings, function(payload_name, artifact_name) {
    relative_path <- paths[[artifact_name]]
    fit_manifest$payloads[[payload_name]]$uri <<- store$object_uri(relative_path)
    store$upload(local_artifacts[[artifact_name]], relative_path)
  })

  fit_manifest$fit_manifest_path <- paths$fit_manifest
  fit_manifest$fit_manifest_uri <- store$object_uri(paths$fit_manifest)
  local_manifest <- tempfile("fit-manifest-", fileext = ".json")
  on.exit(unlink(local_manifest), add = TRUE)
  jsonlite::write_json(
    fit_manifest,
    local_manifest,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null",
    digits = NA
  )
  store$upload(local_manifest, paths$fit_manifest)
  fit_manifest
}

window_manifest_phenotypes <- function(ordered_manifest) {
  required <- c(
    "processing_index",
    "phenotype_id",
    "phenotype_key",
    "modality",
    "p_value"
  )
  missing <- setdiff(required, names(ordered_manifest))
  if (length(missing) > 0L) {
    stop(
      "Ordered manifest is missing checkpoint columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  ordered_manifest |>
    dplyr::select(dplyr::all_of(required)) |>
    purrr::pmap(function(processing_index, phenotype_id, phenotype_key, modality, p_value) {
      list(
        processing_index = checkpoint_scalar_index(
          processing_index,
          "Ordered manifest processing index"
        ),
        phenotype_id = as.character(phenotype_id),
        phenotype_key = as.character(phenotype_key),
        modality = as.character(modality),
        p_value = as.numeric(p_value)
      )
    })
}

new_window_run_manifest <- function(
    analysis_id,
    window_id,
    ordered_manifest,
    input_hashes,
    settings,
    source_hashes = list(unspecified = "unspecified"),
    runtime = list(
      runner_image = "unspecified-runtime",
      base_image_digest = "unspecified-base"
    )
) {
  list(
    schema_version = checkpoint_schema_version,
    analysis_id = analysis_id,
    window_id = window_id,
    phenotypes = window_manifest_phenotypes(ordered_manifest),
    input_hashes = as.list(input_hashes),
    settings = settings,
    source_hashes = as.list(source_hashes),
    runtime = runtime,
    last_committed_index = -1L,
    committed = list(),
    status = "RUNNING",
    failure = NULL,
    failure_history = list(),
    recovery_history = list(),
    attempt = 0L,
    final_outputs = NULL
  )
}

checkpoint_committed_record <- function(fit_manifest) {
  fit_manifest_path <- fit_manifest$fit_manifest_path
  if (is.null(fit_manifest_path)) {
    fit_manifest_path <- phenotype_checkpoint_paths(
      fit_manifest$analysis_id,
      fit_manifest$window_id,
      fit_manifest$phenotype_key
    )$fit_manifest
  }
  list(
    processing_index = checkpoint_scalar_index(
      fit_manifest$processing_index,
      "Checkpoint processing index"
    ),
    phenotype_id = fit_manifest$phenotype_id,
    phenotype_key = fit_manifest$phenotype_key,
    modality = fit_manifest$modality,
    p_value = as.numeric(fit_manifest$p_value),
    status = fit_manifest$status,
    fit_manifest_path = fit_manifest_path
  )
}

advance_window_run_manifest <- function(window_manifest, fit_manifest) {
  processing_index <- checkpoint_scalar_index(
    fit_manifest$processing_index,
    "Checkpoint processing index"
  )
  expected_index <- checkpoint_scalar_index(
    window_manifest$last_committed_index,
    "Window last committed index"
  ) + 1L
  if (!identical(processing_index, expected_index)) {
    stop("Checkpoint index is not the next window index.", call. = FALSE)
  }
  if (!identical(fit_manifest$analysis_id, window_manifest$analysis_id) ||
      !identical(fit_manifest$window_id, window_manifest$window_id)) {
    stop("Checkpoint identity does not match the window manifest.", call. = FALSE)
  }
  if (!fit_manifest$status %in% checkpoint_terminal_statuses) {
    stop("Checkpoint status is not terminal.", call. = FALSE)
  }
  p_value <- fit_manifest$p_value
  if (!is.numeric(p_value) || length(p_value) != 1L || is.na(p_value) ||
      !is.finite(p_value) || p_value < 0 || p_value > 1) {
    stop("Checkpoint p_value must be finite and between zero and one.", call. = FALSE)
  }

  committed_record <- checkpoint_committed_record(fit_manifest)
  window_manifest$committed <- c(window_manifest$committed, list(committed_record))
  window_manifest$last_committed_index <- processing_index
  window_manifest$status <- "RUNNING"
  window_manifest
}

fail_window_run_manifest <- function(window_manifest, condition) {
  if (!inherits(condition, "condition")) {
    stop("Window failure must be an R condition.", call. = FALSE)
  }
  if (!is.null(window_manifest$failure)) {
    window_manifest$failure_history <- c(
      window_manifest$failure_history,
      list(window_manifest$failure)
    )
  }
  next_index <- checkpoint_scalar_index(
    window_manifest$last_committed_index,
    "Window last committed index"
  ) + 1L
  failing_phenotype <- if (next_index >= 0L &&
      next_index < length(window_manifest$phenotypes)) {
    window_manifest$phenotypes[[next_index + 1L]]
  } else {
    NULL
  }
  failure_index <- if (is.null(failing_phenotype)) NULL else next_index
  failure_phenotype_id <- if (is.null(failing_phenotype) ||
      !is.list(failing_phenotype) ||
      !is.character(failing_phenotype$phenotype_id) ||
      length(failing_phenotype$phenotype_id) != 1L) {
    NULL
  } else {
    failing_phenotype$phenotype_id
  }
  window_manifest$failure <- list(
    processing_index = failure_index,
    phenotype_id = failure_phenotype_id,
    class = class(condition)[[1L]],
    message = conditionMessage(condition)
  )
  window_manifest$status <- "FAILED"
  window_manifest
}

checkpoint_order_context <- function(ordered_manifest) {
  required <- c(
    "analysis_id",
    "window_id",
    "phenotype_id",
    "phenotype_key",
    "modality",
    "processing_index"
  )
  missing <- setdiff(required, names(ordered_manifest))
  if (length(missing) > 0L) {
    stop(
      "Ordered manifest is missing resume columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(ordered_manifest) == 0L) {
    stop("Ordered manifest must contain at least one phenotype.", call. = FALSE)
  }

  analysis_ids <- unique(as.character(ordered_manifest$analysis_id))
  window_ids <- unique(as.character(ordered_manifest$window_id))
  if (length(analysis_ids) != 1L || is.na(analysis_ids) || !nzchar(analysis_ids)) {
    stop("Ordered manifest must contain one nonempty analysis_id.", call. = FALSE)
  }
  if (length(window_ids) != 1L || is.na(window_ids) || !nzchar(window_ids)) {
    stop("Ordered manifest must contain one nonempty window_id.", call. = FALSE)
  }
  indexes <- tryCatch(
    purrr::map_int(
      ordered_manifest$processing_index,
      ~ checkpoint_scalar_index(.x, "Ordered manifest processing index")
    ),
    error = function(condition) integer()
  )
  if (!identical(indexes, seq_len(nrow(ordered_manifest)) - 1L)) {
    stop("Ordered manifest processing indexes must be gap-free and zero-based.", call. = FALSE)
  }
  if (anyNA(ordered_manifest$phenotype_key) ||
      any(!nzchar(as.character(ordered_manifest$phenotype_key))) ||
      anyDuplicated(ordered_manifest$phenotype_key)) {
    stop("Ordered manifest phenotype keys must be nonempty and unique.", call. = FALSE)
  }

  list(analysis_id = analysis_ids[[1L]], window_id = window_ids[[1L]])
}

checkpoint_fixed_manifest_path <- function(record) {
  phenotype_checkpoint_paths(
    record$analysis_id[[1L]],
    record$window_id[[1L]],
    record$phenotype_key[[1L]]
  )$fit_manifest
}

checkpoint_canonical_named_json <- function(value) {
  if (!is.list(value)) {
    value <- as.list(value)
  }
  value_names <- names(value)
  if (length(value) == 0L || is.null(value_names) || anyNA(value_names) ||
      any(!nzchar(value_names)) || anyDuplicated(value_names)) {
    return(NULL)
  }
  jsonlite::toJSON(
    value[order(value_names)],
    auto_unbox = TRUE,
    null = "null",
    digits = NA
  )
}

checkpoint_named_object_matches <- function(value, expected) {
  if (is.null(expected)) {
    return(TRUE)
  }
  value_json <- checkpoint_canonical_named_json(value)
  expected_json <- checkpoint_canonical_named_json(expected)
  !is.null(value_json) && !is.null(expected_json) &&
    identical(value_json, expected_json)
}

read_valid_fit_manifest <- function(
    store,
    record,
    expected_input_hashes = NULL,
    expected_settings = NULL,
    expected_source_hashes = NULL,
    expected_runtime = NULL
) {
  fixed_manifest_path <- checkpoint_fixed_manifest_path(record)
  local_manifest <- tempfile("resume-fit-manifest-", fileext = ".json")
  on.exit(unlink(local_manifest), add = TRUE)
  if (!store$object_exists(fixed_manifest_path)) {
    stop("Checkpoint manifest is absent: ", fixed_manifest_path, call. = FALSE)
  }
  store$download(fixed_manifest_path, local_manifest)
  fit_manifest <- jsonlite::read_json(local_manifest, simplifyVector = FALSE)
  validate_fit_manifest_structure(fit_manifest, record)
  if (!checkpoint_named_object_matches(
    fit_manifest$input_hashes,
    expected_input_hashes
  )) {
    stop("Checkpoint input hashes do not match the current analysis.", call. = FALSE)
  }
  if (!checkpoint_named_object_matches(fit_manifest$settings, expected_settings)) {
    stop("Checkpoint settings do not match the current analysis.", call. = FALSE)
  }
  if (!checkpoint_named_object_matches(
    fit_manifest$source_hashes,
    expected_source_hashes
  )) {
    stop("Checkpoint source hashes do not match the current analysis.", call. = FALSE)
  }
  if (!checkpoint_named_object_matches(fit_manifest$runtime, expected_runtime)) {
    stop("Checkpoint runtime identity does not match the current analysis.", call. = FALSE)
  }
  fit_manifest$fit_manifest_path <- fixed_manifest_path
  fit_manifest
}

read_valid_boundary_manifest <- function(
    store,
    record,
    expected_input_hashes = NULL,
    expected_settings = NULL,
    expected_source_hashes = NULL,
    expected_runtime = NULL
) {
  fit_manifest <- read_valid_fit_manifest(
    store,
    record,
    expected_input_hashes,
    expected_settings,
    expected_source_hashes,
    expected_runtime
  )

  if (!identical(fit_manifest$status, "SKIPPED")) {
    payload <- fit_manifest$payloads$susie_fit
    expected_rds_path <- phenotype_checkpoint_paths(
      fit_manifest$analysis_id,
      fit_manifest$window_id,
      fit_manifest$phenotype_key,
      payload$sha256
    )$fit_rds
    if (!identical(payload$path, expected_rds_path)) {
      stop("Checkpoint RDS path does not match its checksum.", call. = FALSE)
    }

    local_rds <- tempfile("resume-fit-", fileext = ".rds")
    on.exit(unlink(local_rds), add = TRUE)
    if (!store$object_exists(payload$path)) {
      stop("Checkpoint RDS payload is absent.", call. = FALSE)
    }
    store$download(payload$path, local_rds)
    if (!identical(sha256_file(local_rds), payload$sha256)) {
      stop("Checkpoint RDS checksum does not match its manifest.", call. = FALSE)
    }
    if (!identical(as.numeric(unname(file.info(local_rds)$size)), as.numeric(payload$bytes))) {
      stop("Checkpoint RDS byte size does not match its manifest.", call. = FALSE)
    }
    fit <- readRDS(local_rds)
    validate_checkpointed_susie_fit(fit)
    if (!identical(fit$converged, fit_manifest$converged)) {
      stop("Checkpoint fit convergence does not match its manifest.", call. = FALSE)
    }
  }

  fit_manifest
}

try_read_valid_boundary_manifest <- function(
    store,
    record,
    expected_input_hashes = NULL,
    expected_settings = NULL,
    expected_source_hashes = NULL,
    expected_runtime = NULL
) {
  tryCatch(
    list(
      valid = TRUE,
      manifest = read_valid_boundary_manifest(
        store,
        record,
        expected_input_hashes,
        expected_settings,
        expected_source_hashes,
        expected_runtime
      ),
      error = NULL
    ),
    error = function(condition) {
      if (inherits(condition, "checkpoint_store_error")) {
        stop(condition)
      }
      list(valid = FALSE, manifest = NULL, error = conditionMessage(condition))
    }
  )
}

checkpoint_cursor_character_matches <- function(value, expected) {
  is.character(value) && length(value) == 1L && !is.na(value) &&
    identical(value, as.character(expected))
}

checkpoint_p_value_matches <- function(value, expected) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || !is.numeric(expected) || length(expected) != 1L ||
      is.na(expected) || !is.finite(expected)) {
    return(FALSE)
  }
  value <- as.numeric(value)
  expected <- as.numeric(expected)
  if (identical(expected, 0)) {
    return(identical(value, 0))
  }
  abs(value - expected) <= abs(expected) * 8 * .Machine$double.eps
}

window_manifest_phenotypes_are_valid <- function(phenotypes, ordered_manifest) {
  if (!is.list(phenotypes) || length(phenotypes) != nrow(ordered_manifest)) {
    return(FALSE)
  }

  for (list_index in seq_len(nrow(ordered_manifest))) {
    record <- phenotypes[[list_index]]
    expected <- ordered_manifest[list_index, , drop = FALSE]
    record_index <- if (is.list(record)) {
      checkpoint_optional_index(record$processing_index)
    } else {
      NULL
    }
    if (is.null(record_index) ||
        !identical(record_index, expected$processing_index[[1L]]) ||
        !checkpoint_cursor_character_matches(
          record$phenotype_key,
          expected$phenotype_key[[1L]]
        ) ||
        !checkpoint_cursor_character_matches(
          record$phenotype_id,
          expected$phenotype_id[[1L]]
        ) ||
        !checkpoint_cursor_character_matches(record$modality, expected$modality[[1L]]) ||
        !checkpoint_p_value_matches(record$p_value, expected$p_value[[1L]])) {
      return(FALSE)
    }
  }
  TRUE
}

window_manifest_commits_are_valid <- function(
    committed,
    ordered_manifest,
    last_committed_index,
    analysis_id,
    window_id,
    recovered_prefix_last_committed_index = -1L
) {
  first_detailed_index <- recovered_prefix_last_committed_index + 1L
  expected_indexes <- if (first_detailed_index <= last_committed_index) {
    seq.int(first_detailed_index, last_committed_index)
  } else {
    integer()
  }
  expected_count <- length(expected_indexes)
  if (!is.list(committed) || length(committed) != expected_count) {
    return(FALSE)
  }
  if (expected_count == 0L) {
    return(TRUE)
  }

  for (list_index in seq_len(expected_count)) {
    record <- committed[[list_index]]
    expected_index <- expected_indexes[[list_index]]
    expected <- ordered_manifest[expected_index + 1L, , drop = FALSE]
    record_index <- if (is.list(record)) {
      checkpoint_optional_index(record$processing_index)
    } else {
      NULL
    }
    if (is.null(record_index) ||
        !identical(record_index, expected$processing_index[[1L]]) ||
        !checkpoint_cursor_character_matches(
          record$phenotype_key,
          expected$phenotype_key[[1L]]
        ) ||
        !checkpoint_cursor_character_matches(
          record$phenotype_id,
          expected$phenotype_id[[1L]]
        ) ||
        !checkpoint_cursor_character_matches(record$modality, expected$modality[[1L]]) ||
        !checkpoint_p_value_matches(record$p_value, expected$p_value[[1L]]) ||
        !is.character(record$status) || length(record$status) != 1L ||
        !record$status %in% checkpoint_terminal_statuses ||
        !is.character(record$fit_manifest_path) ||
        length(record$fit_manifest_path) != 1L ||
        !identical(
          record$fit_manifest_path,
          phenotype_checkpoint_paths(
            analysis_id,
            window_id,
            expected$phenotype_key[[1L]]
          )$fit_manifest
        )) {
      return(FALSE)
    }
  }
  TRUE
}

is_usable_window_manifest <- function(
    window_manifest,
    ordered_manifest,
    context,
    expected_input_hashes = NULL,
    expected_settings = NULL,
    expected_source_hashes = NULL,
    expected_runtime = NULL
) {
  if (!is.list(window_manifest)) {
    return(FALSE)
  }
  schema_version <- checkpoint_optional_index(window_manifest$schema_version)
  last_index <- checkpoint_optional_index(window_manifest$last_committed_index)
  committed <- window_manifest$committed
  if (is.null(committed)) {
    committed <- list()
  }
  status <- window_manifest$status
  valid_status <- is.character(status) && length(status) == 1L &&
    !is.na(status) && status %in% checkpoint_window_statuses
  recovered_prefix <- if (is.null(
    window_manifest$recovered_prefix_last_committed_index
  )) {
    -1L
  } else {
    checkpoint_optional_index(
      window_manifest$recovered_prefix_last_committed_index
    )
  }
  valid_prefix <- !is.null(recovered_prefix) &&
    recovered_prefix >= -1L && !is.null(last_index) &&
    recovered_prefix <= last_index
  inventory_prefix <- if (valid_status && status %in% c("RUNNING", "FAILED")) {
    recovered_prefix
  } else {
    -1L
  }

  !is.null(schema_version) && identical(schema_version, checkpoint_schema_version) &&
    identical(window_manifest$analysis_id, context$analysis_id) &&
    identical(window_manifest$window_id, context$window_id) &&
    valid_status && valid_prefix &&
    checkpoint_named_object_matches(
      window_manifest$input_hashes,
      expected_input_hashes
    ) &&
    checkpoint_named_object_matches(window_manifest$settings, expected_settings) &&
    checkpoint_named_object_matches(
      window_manifest$source_hashes,
      expected_source_hashes
    ) &&
    checkpoint_named_object_matches(window_manifest$runtime, expected_runtime) &&
    !is.null(last_index) &&
    last_index >= -1L && last_index < nrow(ordered_manifest) &&
    window_manifest_phenotypes_are_valid(
      window_manifest$phenotypes,
      ordered_manifest
    ) &&
    window_manifest_commits_are_valid(
      committed,
      ordered_manifest,
      last_index,
      context$analysis_id,
      context$window_id,
      inventory_prefix
    )
}

trim_window_manifest_boundary <- function(window_manifest, invalid_index) {
  committed_indexes <- purrr::map_int(
    window_manifest$committed %||% list(),
    function(record) {
      checkpoint_optional_index(record$processing_index) %||% -1L
    }
  )
  recovered_prefix <- checkpoint_optional_index(
    window_manifest$recovered_prefix_last_committed_index
  ) %||% -1L
  includes_recovered_prefix <- any(committed_indexes <= recovered_prefix)
  if (includes_recovered_prefix) {
    window_manifest$recovered_prefix_last_committed_index <- -1L
  } else if (invalid_index <= recovered_prefix) {
    window_manifest$recovered_prefix_last_committed_index <- invalid_index - 1L
  }
  window_manifest$committed <- purrr::keep(
    window_manifest$committed,
    function(record) {
      record_index <- if (is.list(record)) {
        checkpoint_optional_index(record$processing_index)
      } else {
        NULL
      }
      !is.null(record_index) && record_index < invalid_index
    }
  )
  window_manifest$last_committed_index <- invalid_index - 1L
  window_manifest$status <- "RUNNING"
  window_manifest
}

window_manifest_inventory_prefix <- function(window_manifest) {
  if (!window_manifest$status %in% c("RUNNING", "FAILED") ||
      is.null(window_manifest$recovered_prefix_last_committed_index)) {
    return(-1L)
  }
  checkpoint_scalar_index(
    window_manifest$recovered_prefix_last_committed_index,
    "Window recovered prefix index"
  )
}

window_manifest_cursor_record <- function(window_manifest, processing_index) {
  inventory_prefix <- window_manifest_inventory_prefix(window_manifest)
  list_index <- processing_index - inventory_prefix
  if (list_index < 1L || list_index > length(window_manifest$committed)) {
    return(NULL)
  }
  window_manifest$committed[[list_index]]
}

cursor_record_matches_boundary <- function(cursor_record, boundary_manifest, expected_path) {
  if (!is.list(cursor_record)) {
    return(FALSE)
  }
  cursor_index <- checkpoint_optional_index(cursor_record$processing_index)
  boundary_index <- checkpoint_optional_index(boundary_manifest$processing_index)
  !is.null(cursor_index) && !is.null(boundary_index) &&
    identical(cursor_index, boundary_index) &&
    identical(cursor_record$phenotype_key, boundary_manifest$phenotype_key) &&
    identical(cursor_record$phenotype_id, boundary_manifest$phenotype_id) &&
    identical(cursor_record$modality, boundary_manifest$modality) &&
    checkpoint_p_value_matches(
      cursor_record$p_value,
      boundary_manifest$p_value
    ) &&
    identical(cursor_record$status, boundary_manifest$status) &&
    identical(cursor_record$fit_manifest_path, expected_path)
}

checkpoint_recovery_entry <- function(
    record,
    fit_manifest_path,
    reason,
    attempt = NULL
) {
  concise_reason <- gsub("[[:space:]]+", " ", trimws(as.character(reason)))
  if (nchar(concise_reason) > 500L) {
    concise_reason <- substr(concise_reason, 1L, 500L)
  }
  list(
    processing_index = as.integer(record$processing_index[[1L]]),
    phenotype_id = as.character(record$phenotype_id[[1L]]),
    phenotype_key = as.character(record$phenotype_key[[1L]]),
    fit_manifest_path = fit_manifest_path,
    reason = concise_reason,
    recorded_at = checkpointed_timestamp(),
    attempt = attempt
  )
}

append_window_recovery <- function(window_manifest, entry) {
  window_manifest$recovery_history <- c(
    window_manifest$recovery_history %||% list(),
    list(entry)
  )
  window_manifest
}

resolve_resume_boundary <- function(
    store,
    ordered_manifest,
    window_manifest = NULL,
    expected_input_hashes = NULL,
    expected_settings = NULL,
    expected_source_hashes = NULL,
    expected_runtime = NULL
) {
  context <- checkpoint_order_context(ordered_manifest)
  phenotype_count <- nrow(ordered_manifest)

  if (is_usable_window_manifest(
    window_manifest,
    ordered_manifest,
    context,
    expected_input_hashes,
    expected_settings,
    expected_source_hashes,
    expected_runtime
  )) {
    last_index <- checkpoint_scalar_index(
      window_manifest$last_committed_index,
      "Window last committed index"
    )
    inventory_prefix <- window_manifest_inventory_prefix(window_manifest)
    if (inventory_prefix >= 0L) {
      prefix_record <- ordered_manifest[inventory_prefix + 1L, , drop = FALSE]
      prefix_path <- checkpoint_fixed_manifest_path(prefix_record)
      prefix_boundary <- if (store$object_exists(prefix_path)) {
        try_read_valid_boundary_manifest(
          store,
          prefix_record,
          expected_input_hashes,
          expected_settings,
          expected_source_hashes,
          expected_runtime
        )
      } else {
        list(valid = FALSE, manifest = NULL, error = "Checkpoint manifest is absent.")
      }
      if (!prefix_boundary$valid) {
        attempt <- checkpoint_optional_index(window_manifest$attempt)
        recovery <- checkpoint_recovery_entry(
          prefix_record,
          prefix_path,
          prefix_boundary$error,
          if (is.null(attempt)) NULL else attempt + 1L
        )
        return(list(
          last_committed_index = inventory_prefix - 1L,
          next_index = inventory_prefix,
          recovered = TRUE,
          window_manifest = NULL,
          recovery_history = c(
            window_manifest$recovery_history %||% list(),
            list(recovery)
          )
        ))
      }
    }
    if (last_index >= 0L && last_index != inventory_prefix) {
      record <- ordered_manifest[last_index + 1L, , drop = FALSE]
      expected_path <- checkpoint_fixed_manifest_path(record)
      boundary <- if (store$object_exists(expected_path)) {
        try_read_valid_boundary_manifest(
          store,
          record,
          expected_input_hashes,
          expected_settings,
          expected_source_hashes,
          expected_runtime
        )
      } else {
        list(valid = FALSE, manifest = NULL, error = "Checkpoint manifest is absent.")
      }
      cursor_record <- window_manifest_cursor_record(window_manifest, last_index)
      if (!boundary$valid ||
          !cursor_record_matches_boundary(cursor_record, boundary$manifest, expected_path)) {
        repaired_manifest <- trim_window_manifest_boundary(window_manifest, last_index)
        reason <- if (!boundary$valid) {
          boundary$error
        } else {
          "Window cursor record does not match the boundary checkpoint."
        }
        attempt <- checkpoint_optional_index(window_manifest$attempt)
        repaired_manifest <- append_window_recovery(
          repaired_manifest,
          checkpoint_recovery_entry(
            record,
            expected_path,
            reason,
            if (is.null(attempt)) NULL else attempt + 1L
          )
        )
        return(list(
          last_committed_index = last_index - 1L,
          next_index = last_index,
          recovered = TRUE,
          window_manifest = repaired_manifest
        ))
      }
    }

    next_index <- last_index + 1L
    if (next_index < phenotype_count) {
      next_record <- ordered_manifest[next_index + 1L, , drop = FALSE]
      next_path <- checkpoint_fixed_manifest_path(next_record)
      if (store$object_exists(next_path)) {
        next_boundary <- try_read_valid_boundary_manifest(
          store,
          next_record,
          expected_input_hashes,
          expected_settings,
          expected_source_hashes,
          expected_runtime
        )
        if (next_boundary$valid) {
          advanced_manifest <- advance_window_run_manifest(
            window_manifest,
            next_boundary$manifest
          )
          return(list(
            last_committed_index = next_index,
            next_index = next_index + 1L,
            recovered = TRUE,
            window_manifest = advanced_manifest
          ))
        }
        attempt <- checkpoint_optional_index(window_manifest$attempt)
        audited_manifest <- append_window_recovery(
          window_manifest,
          checkpoint_recovery_entry(
            next_record,
            next_path,
            next_boundary$error,
            if (is.null(attempt)) NULL else attempt + 1L
          )
        )
        return(list(
          last_committed_index = last_index,
          next_index = next_index,
          recovered = TRUE,
          window_manifest = audited_manifest
        ))
      }
    }

    return(list(
      last_committed_index = last_index,
      next_index = next_index,
      recovered = FALSE,
      window_manifest = window_manifest
    ))
  }

  lower_present <- -1L
  upper_absent <- phenotype_count
  while (upper_absent - lower_present > 1L) {
    probe_index <- as.integer(floor((lower_present + upper_absent) / 2))
    probe_record <- ordered_manifest[probe_index + 1L, , drop = FALSE]
    if (store$object_exists(checkpoint_fixed_manifest_path(probe_record))) {
      lower_present <- probe_index
    } else {
      upper_absent <- probe_index
    }
  }

  if (lower_present < 0L) {
    return(list(
      last_committed_index = -1L,
      next_index = 0L,
      recovered = TRUE,
      window_manifest = NULL
    ))
  }

  boundary_record <- ordered_manifest[lower_present + 1L, , drop = FALSE]
  boundary <- try_read_valid_boundary_manifest(
    store,
    boundary_record,
    expected_input_hashes,
    expected_settings,
    expected_source_hashes,
    expected_runtime
  )
  if (!boundary$valid) {
    recovery <- checkpoint_recovery_entry(
      boundary_record,
      checkpoint_fixed_manifest_path(boundary_record),
      boundary$error
    )
    return(list(
      last_committed_index = lower_present - 1L,
      next_index = lower_present,
      recovered = TRUE,
      window_manifest = NULL,
      recovery_history = list(recovery)
    ))
  }

  list(
    last_committed_index = lower_present,
    next_index = lower_present + 1L,
    recovered = TRUE,
    window_manifest = NULL
  )
}

commit_window_outputs <- function(
    store,
    analysis_id,
    window_id,
    local_paths
) {
  required <- c("fit_index", "credible_sets", "lbf_variable", "full_susie")
  if (!is.list(local_paths) || !identical(names(local_paths), required)) {
    stop(
      "Window output paths must contain fit_index and all three Parquet outputs.",
      call. = FALSE
    )
  }
  canonical_paths <- window_checkpoint_paths(analysis_id, window_id)[required]
  records <- purrr::map2(
    local_paths,
    canonical_paths,
    function(local_path, relative_path) {
      record <- checkpointed_payload_record(local_path, relative_path)
      record$uri <- store$object_uri(relative_path)
      record
    }
  )
  purrr::walk2(
    local_paths,
    canonical_paths,
    ~ store$upload(.x, .y)
  )
  records
}

complete_window_run_manifest <- function(
    window_manifest,
    committed_window_outputs
) {
  expected_outputs <- c(
    "fit_index",
    "credible_sets",
    "lbf_variable",
    "full_susie"
  )
  if (!is.list(committed_window_outputs) ||
      !identical(names(committed_window_outputs), expected_outputs)) {
    stop("Committed window outputs are incomplete.", call. = FALSE)
  }
  purrr::iwalk(committed_window_outputs, function(record, name) {
    validate_checkpoint_payload_record(record, paste0("window ", name))
    checkpoint_scalar_character(
      record$uri,
      paste0("Checkpoint window ", name, " URI")
    )
  })
  phenotype_count <- length(window_manifest$phenotypes)
  last_index <- checkpoint_scalar_index(
    window_manifest$last_committed_index,
    "Window last committed index"
  )
  if (!identical(last_index, phenotype_count - 1L) ||
      length(window_manifest$committed) != phenotype_count) {
    stop("Window completion requires every phenotype commit.", call. = FALSE)
  }
  window_manifest$status <- "COMPLETE"
  window_manifest["failure"] <- list(NULL)
  window_manifest$completed_at <- checkpointed_timestamp()
  window_manifest$final_outputs <- committed_window_outputs
  window_manifest
}

write_local_window_run_manifest <- function(window_manifest, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  local_path <- file.path(output_dir, "window_manifest.json")
  jsonlite::write_json(
    window_manifest,
    local_path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null",
    digits = NA
  )
  local_path
}
