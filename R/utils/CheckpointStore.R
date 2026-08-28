# Exact-path checkpoint storage and resume helpers.

checkpoint_schema_version <- 1L

checkpoint_terminal_statuses <- c("CONVERGED", "NONCONVERGED", "SKIPPED")

checkpoint_skip_reasons <- c(
  "NO_USABLE_VARIANTS",
  "ZERO_PHENOTYPE_VARIANCE",
  "NO_ALTERNATIVE_ALLELE",
  "TOO_FEW_ALIGNED_SAMPLES"
)

new_checkpoint_store <- function(root, gsutil = "gsutil") {
  is_gcs <- startsWith(root, "gs://")
  clean_root <- sub("/+$", "", root)
  object_uri <- function(relative_path) paste0(clean_root, "/", relative_path)
  run_gsutil <- function(arguments, operation, uri) {
    status <- system2(gsutil, arguments)
    if (!identical(status, 0L)) {
      stop("GCS ", operation, " failed for: ", uri, call. = FALSE)
    }
    invisible(TRUE)
  }

  list(
    root = clean_root,
    is_gcs = is_gcs,
    object_exists = function(relative_path) {
      if (is_gcs) {
        return(identical(system2(gsutil, c("-q", "stat", object_uri(relative_path))), 0L))
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
        stop("Filesystem upload failed for: ", destination, call. = FALSE)
      }
      invisible(TRUE)
    },
    download = function(relative_path, local_path) {
      dir.create(dirname(local_path), recursive = TRUE, showWarnings = FALSE)
      if (is_gcs) {
        return(run_gsutil(c("-q", "cp", object_uri(relative_path), local_path), "download", object_uri(relative_path)))
      }
      source_path <- file.path(clean_root, relative_path)
      if (!file.copy(source_path, local_path, overwrite = TRUE)) {
        stop("Filesystem download failed for: ", source_path, call. = FALSE)
      }
      invisible(TRUE)
    },
    object_uri = object_uri
  )
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
  numeric_value <- suppressWarnings(as.numeric(value))
  if (length(numeric_value) != 1L || is.na(numeric_value) ||
      !is.finite(numeric_value) || numeric_value != floor(numeric_value)) {
    stop(label, " must be one integer.", call. = FALSE)
  }
  as.integer(numeric_value)
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
    expected_index <- as.integer(record$processing_index[[1L]])
    if (!identical(processing_index, expected_index)) {
      stop("Checkpoint processing index does not match the ordered manifest.", call. = FALSE)
    }
  }

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
    return(invisible(TRUE))
  }

  if (!is.logical(fit_manifest$converged) ||
      length(fit_manifest$converged) != 1L || is.na(fit_manifest$converged)) {
    stop("Checkpoint convergence must be one logical value.", call. = FALSE)
  }
  expected_convergence <- identical(status, "CONVERGED")
  if (!identical(fit_manifest$converged, expected_convergence)) {
    stop("Checkpoint status and convergence do not match.", call. = FALSE)
  }
  required_payloads <- c("susie_fit", "credible_sets", "lbf_variable", "full_susie")
  purrr::walk(
    required_payloads,
    ~ validate_checkpoint_payload_record(fit_manifest$payloads[[.x]], .x)
  )
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
    null = "null"
  )
  store$upload(local_manifest, paths$fit_manifest)
  fit_manifest
}

window_manifest_phenotypes <- function(ordered_manifest) {
  required <- c("processing_index", "phenotype_id", "phenotype_key", "modality")
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
    purrr::pmap(function(processing_index, phenotype_id, phenotype_key, modality) {
      list(
        processing_index = as.integer(processing_index),
        phenotype_id = as.character(phenotype_id),
        phenotype_key = as.character(phenotype_key),
        modality = as.character(modality)
      )
    })
}

new_window_run_manifest <- function(
    analysis_id,
    window_id,
    ordered_manifest,
    input_hashes,
    settings
) {
  list(
    schema_version = checkpoint_schema_version,
    analysis_id = analysis_id,
    window_id = window_id,
    phenotypes = window_manifest_phenotypes(ordered_manifest),
    input_hashes = as.list(input_hashes),
    settings = settings,
    last_committed_index = -1L,
    committed = list(),
    status = "RUNNING",
    failure = NULL,
    failure_history = list(),
    final_outputs = NULL
  )
}

advance_window_run_manifest <- function(window_manifest, fit_manifest) {
  processing_index <- checkpoint_scalar_index(
    fit_manifest$processing_index,
    "Checkpoint processing index"
  )
  expected_index <- as.integer(window_manifest$last_committed_index) + 1L
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

  fit_manifest_path <- fit_manifest$fit_manifest_path
  if (is.null(fit_manifest_path)) {
    fit_manifest_path <- phenotype_checkpoint_paths(
      fit_manifest$analysis_id,
      fit_manifest$window_id,
      fit_manifest$phenotype_key
    )$fit_manifest
  }
  committed_record <- list(
    processing_index = processing_index,
    phenotype_id = fit_manifest$phenotype_id,
    phenotype_key = fit_manifest$phenotype_key,
    modality = fit_manifest$modality,
    status = fit_manifest$status,
    fit_manifest_path = fit_manifest_path
  )
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
  window_manifest$failure <- list(
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
  indexes <- suppressWarnings(as.integer(ordered_manifest$processing_index))
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

read_valid_boundary_manifest <- function(store, record) {
  fixed_manifest_path <- checkpoint_fixed_manifest_path(record)
  local_manifest <- tempfile("resume-fit-manifest-", fileext = ".json")
  on.exit(unlink(local_manifest), add = TRUE)
  store$download(fixed_manifest_path, local_manifest)
  fit_manifest <- jsonlite::read_json(local_manifest, simplifyVector = FALSE)
  validate_fit_manifest_structure(fit_manifest, record)
  fit_manifest$fit_manifest_path <- fixed_manifest_path

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

try_read_valid_boundary_manifest <- function(store, record) {
  tryCatch(
    list(valid = TRUE, manifest = read_valid_boundary_manifest(store, record), error = NULL),
    error = function(condition) {
      list(valid = FALSE, manifest = NULL, error = conditionMessage(condition))
    }
  )
}

window_manifest_keys <- function(window_manifest) {
  phenotypes <- window_manifest$phenotypes
  if (!is.list(phenotypes)) {
    return(character())
  }
  purrr::map_chr(phenotypes, function(record) {
    if (!is.list(record) || is.null(record$phenotype_key)) {
      return(NA_character_)
    }
    as.character(record$phenotype_key)
  })
}

is_usable_window_manifest <- function(window_manifest, ordered_manifest, context) {
  if (!is.list(window_manifest)) {
    return(FALSE)
  }
  schema_version <- suppressWarnings(as.integer(window_manifest$schema_version))
  last_index <- suppressWarnings(as.integer(window_manifest$last_committed_index))
  expected_keys <- as.character(ordered_manifest$phenotype_key)
  committed <- window_manifest$committed
  if (is.null(committed)) {
    committed <- list()
  }

  length(schema_version) == 1L && !is.na(schema_version) &&
    identical(schema_version, checkpoint_schema_version) &&
    identical(window_manifest$analysis_id, context$analysis_id) &&
    identical(window_manifest$window_id, context$window_id) &&
    length(last_index) == 1L && !is.na(last_index) &&
    last_index >= -1L && last_index < nrow(ordered_manifest) &&
    identical(window_manifest_keys(window_manifest), expected_keys) &&
    length(committed) >= last_index + 1L
}

trim_window_manifest_boundary <- function(window_manifest, invalid_index) {
  retained_count <- max(0L, invalid_index)
  window_manifest$committed <- utils::head(window_manifest$committed, retained_count)
  window_manifest$last_committed_index <- invalid_index - 1L
  window_manifest$status <- "RUNNING"
  window_manifest
}

cursor_record_matches_boundary <- function(cursor_record, boundary_manifest, expected_path) {
  if (!is.list(cursor_record)) {
    return(FALSE)
  }
  identical(as.integer(cursor_record$processing_index), as.integer(boundary_manifest$processing_index)) &&
    identical(cursor_record$phenotype_key, boundary_manifest$phenotype_key) &&
    identical(cursor_record$status, boundary_manifest$status) &&
    identical(cursor_record$fit_manifest_path, expected_path)
}

resolve_resume_boundary <- function(store, ordered_manifest, window_manifest = NULL) {
  context <- checkpoint_order_context(ordered_manifest)
  phenotype_count <- nrow(ordered_manifest)

  if (is_usable_window_manifest(window_manifest, ordered_manifest, context)) {
    last_index <- as.integer(window_manifest$last_committed_index)
    if (last_index >= 0L) {
      record <- ordered_manifest[last_index + 1L, , drop = FALSE]
      boundary <- try_read_valid_boundary_manifest(store, record)
      expected_path <- checkpoint_fixed_manifest_path(record)
      cursor_record <- window_manifest$committed[[last_index + 1L]]
      if (!boundary$valid ||
          !cursor_record_matches_boundary(cursor_record, boundary$manifest, expected_path)) {
        repaired_manifest <- trim_window_manifest_boundary(window_manifest, last_index)
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
        next_boundary <- try_read_valid_boundary_manifest(store, next_record)
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
  boundary <- try_read_valid_boundary_manifest(store, boundary_record)
  if (!boundary$valid) {
    return(list(
      last_committed_index = lower_present - 1L,
      next_index = lower_present,
      recovered = TRUE,
      window_manifest = NULL
    ))
  }

  list(
    last_committed_index = lower_present,
    next_index = lower_present + 1L,
    recovered = TRUE,
    window_manifest = NULL
  )
}
