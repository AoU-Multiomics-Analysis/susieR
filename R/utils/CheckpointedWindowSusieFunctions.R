# Canonical manifest and identity helpers for one-window univariate SuSiE.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

checkpointed_susie_settings <- function() {
  list(
    L = 10L,
    estimate_residual_variance = TRUE,
    estimate_prior_variance = TRUE,
    scaled_prior_variance = 0.1,
    compute_univariate_zscore = TRUE,
    min_abs_corr = 0.5
  )
}

checkpointed_qr_tolerance <- 1e-7

checkpoint_phenotype_key <- function(window_id, modality, phenotype_id) {
  digest::digest(
    paste(window_id, modality, phenotype_id, sep = "\n"),
    algo = "sha256",
    serialize = FALSE
  )
}

validate_window_phenotype_manifest <- function(
    manifest,
    window_id,
    phenotype_data = NULL
) {
  required <- c("window_id", "phenotype_id", "modality", "phenotype_file", "p_value")
  missing <- setdiff(required, names(manifest))
  if (length(missing) > 0L) {
    stop(
      "Prepared manifest is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  validated <- manifest |>
    dplyr::transmute(
      window_id = as.character(.data$window_id),
      phenotype_id = as.character(.data$phenotype_id),
      modality = as.character(.data$modality),
      phenotype_file = as.character(.data$phenotype_file),
      p_value = suppressWarnings(as.numeric(.data$p_value))
    )

  identifier_columns <- c("window_id", "phenotype_id", "modality", "phenotype_file")
  has_empty_identifier <- purrr::map_lgl(
    validated[identifier_columns],
    ~ any(is.na(.x) | !nzchar(.x))
  ) |>
    any()
  if (has_empty_identifier) {
    stop("Prepared manifest identifiers cannot be empty.", call. = FALSE)
  }
  if (anyNA(validated$p_value) || any(!is.finite(validated$p_value)) ||
      any(validated$p_value < 0 | validated$p_value > 1)) {
    stop(
      "Prepared manifest p_value values must be finite and between zero and one.",
      call. = FALSE
    )
  }
  if (anyDuplicated(validated$phenotype_id)) {
    stop("Prepared manifest contains duplicate phenotype_id values.", call. = FALSE)
  }
  manifest_windows <- unique(validated$window_id)
  if (length(manifest_windows) != 1L || !identical(manifest_windows, window_id)) {
    stop("Prepared manifest window_id does not match the requested window.", call. = FALSE)
  }
  phenotype_files <- unique(validated$phenotype_file)
  if (length(phenotype_files) != 1L) {
    stop("Prepared manifest must contain one phenotype_file value.", call. = FALSE)
  }
  if (!is.null(phenotype_data) &&
      !identical(basename(phenotype_files), basename(phenotype_data))) {
    stop(
      "Prepared manifest phenotype_file must match the configured phenotype data file.",
      call. = FALSE
    )
  }

  validated |>
    dplyr::arrange(.data$p_value, .data$modality, .data$phenotype_id) |>
    dplyr::mutate(
      processing_index = dplyr::row_number() - 1L,
      phenotype_key = purrr::pmap_chr(
        list(.data$window_id, .data$modality, .data$phenotype_id),
        checkpoint_phenotype_key
      )
    )
}

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

build_checkpoint_analysis_id <- function(
    input_hashes,
    ordered_manifest,
    settings,
    container_digest,
    wrapper_hash,
    source_hashes = NULL,
    runtime_image = "unspecified-runtime-image"
) {
  if (is.null(source_hashes)) {
    source_hashes <- c(runner_wrapper = wrapper_hash)
  }
  canonical_payload <- list(
    input_hashes = as.list(input_hashes[order(names(input_hashes))]),
    phenotypes = ordered_manifest |>
      dplyr::select(
        processing_index,
        window_id,
        phenotype_id,
        modality,
        p_value,
        phenotype_key
      ),
    settings = settings[order(names(settings))],
    runtime = list(
      runner_image = runtime_image,
      base_image_digest = container_digest
    ),
    source_hashes = as.list(source_hashes[order(names(source_hashes))])
  )
  canonical_json <- jsonlite::toJSON(
    canonical_payload,
    auto_unbox = TRUE,
    dataframe = "rows",
    null = "null",
    digits = NA
  )
  digest::digest(canonical_json, algo = "sha256", serialize = FALSE)
}

validate_checkpointed_susie_fit <- function(fit) {
  if (!inherits(fit, "susie")) {
    stop("Checkpoint fit must inherit from class susie.", call. = FALSE)
  }

  required_fields <- c("pip", "alpha", "mu", "mu2", "variant_id", "converged")
  missing_fields <- setdiff(required_fields, names(fit))
  if (length(missing_fields) > 0L) {
    stop(
      "Checkpoint fit is missing required fields: ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  variant_id <- fit$variant_id
  if (!is.character(variant_id) ||
      anyNA(variant_id) || any(!nzchar(variant_id))) {
    stop("Checkpoint variant IDs must be nonempty strings.", call. = FALSE)
  }
  if (anyDuplicated(variant_id)) {
    stop("Checkpoint variant IDs must be unique.", call. = FALSE)
  }

  pip <- fit$pip
  if (!is.numeric(pip) || !is.null(dim(pip)) || length(pip) != length(variant_id) ||
      anyNA(pip) || any(!is.finite(pip)) || any(pip < 0 | pip > 1)) {
    stop(
      "Checkpoint PIP values must be a finite numeric vector between zero and one.",
      call. = FALSE
    )
  }

  model_arrays <- fit[c("alpha", "mu", "mu2")]
  valid_arrays <- purrr::map_lgl(
    model_arrays,
    ~ is.matrix(.x) && is.numeric(.x) && all(is.finite(.x))
  )
  if (!all(valid_arrays)) {
    stop("Checkpoint model arrays must contain finite numeric values.", call. = FALSE)
  }
  array_dimensions <- purrr::map(model_arrays, dim)
  positive_dimensions <- length(array_dimensions[[1L]]) == 2L &&
    all(array_dimensions[[1L]] > 0L)
  if (!positive_dimensions) {
    stop("Checkpoint model arrays must have positive dimensions.", call. = FALSE)
  }
  if (!all(purrr::map_lgl(array_dimensions, identical, array_dimensions[[1L]])) ||
      ncol(model_arrays[[1L]]) != length(variant_id)) {
    stop(
      "Checkpoint alpha, mu, and mu2 arrays must have the same dimensions as the variant data.",
      call. = FALSE
    )
  }

  if (!is.logical(fit$converged) || length(fit$converged) != 1L ||
      is.na(fit$converged)) {
    stop("Checkpoint convergence must be one logical value.", call. = FALSE)
  }

  TRUE
}

checkpointed_numeric_matrix <- function(data, columns, label) {
  numeric_columns <- purrr::map_lgl(data[columns], is.numeric)
  if (!all(numeric_columns)) {
    stop(label, " sample values must be numeric.", call. = FALSE)
  }
  matrix <- as.matrix(data[columns])
  storage.mode(matrix) <- "double"
  matrix
}

read_checkpointed_window_dosage <- function(path, impute = TRUE) {
  dosage <- data.table::fread(
    path,
    check.names = FALSE,
    data.table = FALSE
  )
  chromosome_column <- intersect(c("CHROM", "#CHROM"), names(dosage))
  if (length(chromosome_column) != 1L) {
    stop("Prepared dosage must contain one CHROM or #CHROM column.", call. = FALSE)
  }
  required <- c("POS", "REF", "ALT")
  missing <- setdiff(required, names(dosage))
  if (length(missing) > 0L) {
    stop(
      "Prepared dosage is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  metadata_columns <- c(chromosome_column, required)
  sample_ids <- setdiff(names(dosage), metadata_columns)
  if (length(sample_ids) == 0L || any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
    stop("Prepared dosage sample IDs must be unique nonempty strings.", call. = FALSE)
  }
  genotype <- checkpointed_numeric_matrix(dosage, sample_ids, "Prepared dosage")
  genotype[genotype == -1 & !is.na(genotype)] <- NA_real_
  if (any(!is.finite(genotype) & !is.na(genotype))) {
    stop("Prepared dosage contains nonfinite values.", call. = FALSE)
  }
  invalid_dosage <- is.finite(genotype) & (genotype < 0 | genotype > 2)
  if (any(invalid_dosage)) {
    stop("Prepared dosage values must be between zero and two.", call. = FALSE)
  }

  if (isTRUE(impute)) {
    imputation_means <- rowMeans(genotype, na.rm = TRUE)
    missing_cells <- which(is.na(genotype), arr.ind = TRUE)
    if (nrow(missing_cells) > 0L) {
      replacement_values <- imputation_means[missing_cells[, "row"]]
      finite_replacements <- is.finite(replacement_values)
      genotype[missing_cells[finite_replacements, , drop = FALSE]] <-
        replacement_values[finite_replacements]
    }
  }

  variant_info <- dosage |>
    dplyr::transmute(
      CHROM = as.character(.data[[chromosome_column]]),
      POS = as.integer(.data$POS),
      REF = as.character(.data$REF),
      ALT = as.character(.data$ALT),
      variant_id = paste(.data$CHROM, .data$POS, .data$REF, .data$ALT, sep = "_")
    )
  if (anyNA(variant_info) || any(!nzchar(variant_info$variant_id)) ||
      anyDuplicated(variant_info$variant_id)) {
    stop("Prepared dosage variant metadata must be complete and unique.", call. = FALSE)
  }

  rownames(genotype) <- variant_info$variant_id
  colnames(genotype) <- sample_ids
  list(genotype = genotype, variant_info = variant_info)
}

checkpointed_column_name <- function(names, choices, label) {
  matches <- intersect(choices, names)
  if (length(matches) != 1L) {
    stop("Prepared phenotype data must contain one ", label, " column.", call. = FALSE)
  }
  matches[[1L]]
}

read_checkpointed_window_phenotypes <- function(path, ordered_manifest) {
  phenotype_data <- data.table::fread(
    path,
    check.names = FALSE,
    data.table = FALSE
  )
  chromosome_column <- checkpointed_column_name(
    names(phenotype_data),
    c("chromosome", "chrom", "CHROM", "#CHROM", "chr", "#chr"),
    "chromosome"
  )
  start_column <- checkpointed_column_name(
    names(phenotype_data),
    c("start", "START"),
    "start"
  )
  end_column <- checkpointed_column_name(
    names(phenotype_data),
    c("end", "END"),
    "end"
  )
  phenotype_id_column <- checkpointed_column_name(
    names(phenotype_data),
    c("phenotype_id"),
    "phenotype_id"
  )
  metadata_columns <- c(
    chromosome_column,
    start_column,
    end_column,
    phenotype_id_column
  )
  sample_ids <- setdiff(names(phenotype_data), metadata_columns)
  if (length(sample_ids) == 0L || any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
    stop("Prepared phenotype sample IDs must be unique nonempty strings.", call. = FALSE)
  }
  if (anyDuplicated(phenotype_data[[phenotype_id_column]])) {
    stop("Prepared phenotype data contains duplicate phenotype IDs.", call. = FALSE)
  }

  requested_ids <- ordered_manifest$phenotype_id
  selected_rows <- match(requested_ids, phenotype_data[[phenotype_id_column]])
  if (anyNA(selected_rows)) {
    stop(
      "Prepared phenotype data is missing manifest phenotypes: ",
      paste(requested_ids[is.na(selected_rows)], collapse = ", "),
      call. = FALSE
    )
  }
  selected <- phenotype_data[selected_rows, , drop = FALSE]
  value_matrix <- checkpointed_numeric_matrix(selected, sample_ids, "Prepared phenotype")
  values <- purrr::map(seq_len(nrow(value_matrix)), function(index) {
    stats::setNames(as.numeric(value_matrix[index, ]), sample_ids)
  }) |>
    stats::setNames(requested_ids)
  metadata <- selected |>
    dplyr::transmute(
      chromosome = as.character(.data[[chromosome_column]]),
      start = as.integer(.data[[start_column]]),
      end = as.integer(.data[[end_column]]),
      phenotype_id = as.character(.data[[phenotype_id_column]])
    )
  list(values = values, metadata = metadata)
}

read_checkpointed_covariates <- function(paths, labels) {
  paths <- as.character(paths)
  labels <- as.character(labels)
  if (length(paths) == 0L || length(paths) != length(labels)) {
    stop("Covariate paths and labels must have the same nonzero length.", call. = FALSE)
  }
  if (anyNA(labels) || any(!nzchar(labels))) {
    stop("Covariate labels must be nonempty strings.", call. = FALSE)
  }

  covariates <- purrr::map2(paths, labels, function(path, label) {
    data <- data.table::fread(
      path,
      header = TRUE,
      check.names = FALSE,
      data.table = FALSE
    )
    if (ncol(data) < 2L) {
      stop("Each covariate file must contain IDs and sample values.", call. = FALSE)
    }
    covariate_id_column <- names(data)[[1L]]
    sample_ids <- names(data)[-1L]
    covariate_ids <- as.character(data[[covariate_id_column]])
    if (anyNA(covariate_ids) || any(!nzchar(covariate_ids)) ||
        anyDuplicated(covariate_ids)) {
      stop("Covariate IDs within each file must be unique nonempty strings.", call. = FALSE)
    }
    if (any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
      stop("Covariate sample IDs must be unique nonempty strings.", call. = FALSE)
    }
    values <- checkpointed_numeric_matrix(data, sample_ids, "Covariate")
    if (any(!is.finite(values))) {
      stop("Covariate values must be finite.", call. = FALSE)
    }
    rownames(values) <- covariate_ids
    colnames(values) <- sample_ids
    list(
      path = path,
      label = label,
      covariate_ids = covariate_ids,
      sample_ids = sample_ids,
      values = values
    )
  })
  stats::setNames(covariates, sprintf("covariate_%03d", seq_along(covariates)))
}

checkpointed_reduce_covariate_matrix <- function(sample_by_covariate) {
  sample_by_covariate <- as.matrix(sample_by_covariate)
  storage.mode(sample_by_covariate) <- "double"
  covariate_ids <- colnames(sample_by_covariate)
  if (is.null(covariate_ids)) {
    covariate_ids <- sprintf("covariate_%03d", seq_len(ncol(sample_by_covariate)))
    colnames(sample_by_covariate) <- covariate_ids
  }
  if (any(!is.finite(sample_by_covariate))) {
    stop("Applicable covariate design must be finite.", call. = FALSE)
  }

  variances <- if (ncol(sample_by_covariate) == 0L) {
    numeric()
  } else {
    apply(sample_by_covariate, 2L, stats::var)
  }
  nonconstant <- is.finite(variances) & variances > 0
  candidate <- sample_by_covariate[, nonconstant, drop = FALSE]
  retained <- integer()
  if (ncol(candidate) > 0L) {
    for (column_index in seq_len(ncol(candidate))) {
      trial <- candidate[, c(retained, column_index), drop = FALSE]
      trial_design <- cbind(intercept = 1, trial)
      trial_rank <- qr(
        trial_design,
        tol = checkpointed_qr_tolerance,
        LAPACK = FALSE
      )$rank
      if (trial_rank > 1L + length(retained)) {
        retained <- c(retained, column_index)
      }
    }
  }
  full_rank <- candidate[, retained, drop = FALSE]
  if (ncol(full_rank) == 0L) {
    full_rank <- matrix(
      numeric(nrow(sample_by_covariate) * 0L),
      nrow = nrow(sample_by_covariate),
      dimnames = list(rownames(sample_by_covariate), character())
    )
  }
  attr(full_rank, "removed_covariates") <- covariate_ids[
    !covariate_ids %in% colnames(full_rank)
  ]
  attr(full_rank, "input_covariates") <- covariate_ids
  full_rank
}

applicable_covariate_matrix <- function(covariates, modality, sample_ids) {
  modality <- as.character(modality)
  sample_ids <- as.character(sample_ids)
  if (length(modality) != 1L || is.na(modality) || !nzchar(modality)) {
    stop("Phenotype modality must be one nonempty string.", call. = FALSE)
  }
  if (anyNA(sample_ids) || any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
    stop("Aligned sample IDs must be unique nonempty strings.", call. = FALSE)
  }
  selected <- purrr::keep(
    covariates,
    ~ .x$label %in% c("shared", modality)
  )
  if (length(selected) == 0L) {
    stop("Phenotype modality has no applicable covariate input.", call. = FALSE)
  }
  purrr::walk(selected, function(covariate) {
    missing_samples <- setdiff(sample_ids, covariate$sample_ids)
    if (length(missing_samples) > 0L) {
      stop("Applicable covariates are missing aligned samples.", call. = FALSE)
    }
  })

  combined <- purrr::map(selected, function(covariate) {
    covariate$values[, sample_ids, drop = FALSE]
  }) |>
    do.call(what = rbind)
  covariate_ids <- rownames(combined)
  if (anyDuplicated(covariate_ids)) {
    stop("Applicable files contain duplicate covariate IDs.", call. = FALSE)
  }
  sample_by_covariate <- t(combined)
  checkpointed_reduce_covariate_matrix(sample_by_covariate)
}

read_checkpointed_sample_allowlist <- function(path) {
  allowlist <- data.table::fread(
    path,
    check.names = FALSE,
    data.table = FALSE
  )
  if (ncol(allowlist) != 1L || is.null(names(allowlist)) ||
      is.na(names(allowlist)[[1L]]) || !nzchar(names(allowlist)[[1L]])) {
    stop("Sample allowlist must contain exactly one named sample-ID column.", call. = FALSE)
  }
  sample_ids <- as.character(allowlist[[1L]])
  if (anyNA(sample_ids) || any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
    stop("Sample allowlist IDs must be unique nonempty strings.", call. = FALSE)
  }
  sample_ids
}

load_checkpointed_window_inputs <- function(config, ordered_manifest) {
  required <- c(
    "window_dosage",
    "phenotype_data",
    "covariate_files",
    "covariate_modalities"
  )
  missing <- setdiff(required, names(config))
  if (length(missing) > 0L) {
    stop(
      "Model input configuration is missing fields: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  dosage <- read_checkpointed_window_dosage(config$window_dosage, impute = FALSE)
  phenotype_data <- read_checkpointed_window_phenotypes(
    config$phenotype_data,
    ordered_manifest
  )
  covariates <- read_checkpointed_covariates(
    config$covariate_files,
    config$covariate_modalities
  )

  dosage_samples <- colnames(dosage$genotype)
  manifest_modalities <- unique(ordered_manifest$modality)
  applicable_covariates <- purrr::keep(
    covariates,
    ~ .x$label %in% c("shared", manifest_modalities)
  )
  purrr::walk(manifest_modalities, function(modality) {
    if (!any(purrr::map_lgl(
      covariates,
      ~ .x$label %in% c("shared", modality)
    ))) {
      stop("Phenotype modality has no applicable covariate input.", call. = FALSE)
    }
  })
  available_sample_sets <- c(
    list(names(phenotype_data$values[[1L]])),
    purrr::map(applicable_covariates, "sample_ids")
  )
  keep_samples <- config$keep_samples
  has_allowlist <- !is.null(keep_samples) && length(keep_samples) == 1L &&
    !is.na(keep_samples) && nzchar(keep_samples)
  if (has_allowlist) {
    available_sample_sets <- c(
      available_sample_sets,
      list(read_checkpointed_sample_allowlist(keep_samples))
    )
  }
  aligned_samples <- dosage_samples[
    purrr::map_lgl(dosage_samples, function(sample_id) {
      all(purrr::map_lgl(available_sample_sets, ~ sample_id %in% .x))
    })
  ]
  genotype <- dosage$genotype[, aligned_samples, drop = FALSE]
  phenotypes <- purrr::map(
    phenotype_data$values,
    ~ .x[aligned_samples]
  )
  list(
    genotype = genotype,
    variant_info = dosage$variant_info,
    phenotypes = phenotypes,
    phenotype_metadata = phenotype_data$metadata,
    covariates = covariates,
    sample_ids = aligned_samples,
    raw_dosage_samples = as.integer(length(dosage_samples)),
    overlap_samples = as.integer(length(aligned_samples))
  )
}

rank_inverse_normal <- function(values) {
  stats::qnorm(
    (rank(values, ties.method = "average") - 0.5) / length(values)
  )
}

checkpointed_has_numerical_variation <- function(values, reference_values) {
  values <- as.numeric(values)
  reference_values <- as.numeric(reference_values)
  if (length(values) < 2L || length(values) != length(reference_values) ||
      any(!is.finite(values)) || any(!is.finite(reference_values))) {
    return(FALSE)
  }
  centered_values <- values - mean(values)
  centered_reference <- reference_values - mean(reference_values)
  residual_norm <- sqrt(sum(centered_values^2))
  reference_norm <- sqrt(sum(centered_reference^2))
  numerical_error_bound <-
    100 * length(values) * .Machine$double.eps * reference_norm
  is.finite(residual_norm) && is.finite(reference_norm) &&
    reference_norm > 0 && residual_norm > numerical_error_bound
}

checkpointed_skipped_fit <- function(reason, qc) {
  list(status = "SKIPPED", fit = NULL, reason = reason, qc = qc)
}

checkpointed_matrix_identity_payload <- function(matrix) {
  matrix <- as.matrix(matrix)
  list(
    dimensions = as.integer(dim(matrix)),
    column_names = colnames(matrix) %||% character(),
    values_by_row = as.numeric(t(matrix))
  )
}

checkpointed_sample_mask_id <- function(sample_ids) {
  digest::digest(
    paste(as.character(sample_ids), collapse = "\n"),
    algo = "sha256",
    serialize = FALSE
  )
}

checkpointed_design_identity <- function(sample_ids, design) {
  payload <- list(
    ordered_sample_ids = as.character(sample_ids),
    design = checkpointed_matrix_identity_payload(design)
  )
  digest::digest(
    jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null", digits = NA),
    algo = "sha256",
    serialize = FALSE
  )
}

checkpointed_genotype_cache_key <- function(
    sample_ids,
    retained_variant_ids,
    design
) {
  payload <- list(
    ordered_sample_ids = as.character(sample_ids),
    retained_variant_ids = as.character(retained_variant_ids),
    retained_covariate_columns = colnames(design) %||% character(),
    design = checkpointed_matrix_identity_payload(design)
  )
  digest::digest(
    jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null", digits = NA),
    algo = "sha256",
    serialize = FALSE
  )
}

build_checkpointed_window_cache_plan <- function(model_inputs, ordered_manifest) {
  required_inputs <- c(
    "genotype",
    "variant_info",
    "phenotypes",
    "covariates",
    "sample_ids"
  )
  missing_inputs <- setdiff(required_inputs, names(model_inputs))
  if (length(missing_inputs) > 0L) {
    stop(
      "Cache planning inputs are missing fields: ",
      paste(missing_inputs, collapse = ", "),
      call. = FALSE
    )
  }

  phenotype_plans <- purrr::map(
    seq_len(nrow(ordered_manifest)),
    function(row_index) {
      record <- ordered_manifest[row_index, , drop = FALSE]
      phenotype_id <- record$phenotype_id[[1L]]
      phenotype <- model_inputs$phenotypes[[phenotype_id]]
      if (is.null(phenotype)) {
        stop("Cache planning phenotype is absent: ", phenotype_id, call. = FALSE)
      }
      aligned_values <- phenotype[model_inputs$sample_ids]
      retained_sample_ids <- model_inputs$sample_ids[is.finite(aligned_values)]
      covariates <- applicable_covariate_matrix(
        model_inputs$covariates,
        record$modality[[1L]],
        retained_sample_ids
      )
      design <- cbind(intercept = 1, covariates)
      rownames(design) <- retained_sample_ids
      design_id <- checkpointed_design_identity(retained_sample_ids, design)
      sample_mask_id <- checkpointed_sample_mask_id(retained_sample_ids)
      preparation_key <- digest::digest(
        paste(design_id, sample_mask_id, sep = "\n"),
        algo = "sha256",
        serialize = FALSE
      )
      list(
        processing_index = as.integer(record$processing_index[[1L]]),
        phenotype_id = phenotype_id,
        modality = record$modality[[1L]],
        retained_sample_ids = retained_sample_ids,
        covariates = covariates,
        design = design,
        design_id = design_id,
        sample_mask_id = sample_mask_id,
        preparation_key = preparation_key,
        removed_covariates = attr(covariates, "removed_covariates") %||%
          character()
      )
    }
  ) |>
    stats::setNames(ordered_manifest$phenotype_id)
  list(phenotypes = phenotype_plans)
}

checkpointed_preparation_qc <- function(
    raw_dosage_samples,
    overlap_samples,
    sample_ids,
    input_variant_ids,
    retained_variant_ids,
    plan,
    cache_key
) {
  list(
    input_samples = as.integer(overlap_samples),
    retained_samples = as.integer(length(sample_ids)),
    raw_dosage_samples = as.integer(raw_dosage_samples),
    overlap_samples = as.integer(overlap_samples),
    phenotype_retained_samples = as.integer(length(sample_ids)),
    input_variants = as.integer(length(input_variant_ids)),
    retained_variants = as.integer(length(retained_variant_ids)),
    removed_covariates = as.character(plan$removed_covariates),
    removed_variant_ids = input_variant_ids[
      !input_variant_ids %in% retained_variant_ids
    ],
    design_id = plan$design_id,
    sample_mask_id = plan$sample_mask_id,
    cache_key = cache_key
  )
}

prepare_checkpointed_window_genotype <- function(
    genotype,
    variant_ids,
    plan,
    raw_dosage_samples,
    overlap_samples
) {
  genotype <- as.matrix(genotype)
  storage.mode(genotype) <- "double"
  variant_ids <- as.character(variant_ids)
  sample_ids <- as.character(plan$retained_sample_ids)
  if (nrow(genotype) != length(variant_ids) || is.null(colnames(genotype)) ||
      any(!sample_ids %in% colnames(genotype))) {
    stop("Genotype and cache-plan dimensions must align.", call. = FALSE)
  }
  if (anyNA(variant_ids) || any(!nzchar(variant_ids)) || anyDuplicated(variant_ids)) {
    stop("Model variant IDs must be unique nonempty strings.", call. = FALSE)
  }
  genotype <- genotype[, sample_ids, drop = FALSE]
  design <- as.matrix(plan$design)
  storage.mode(design) <- "double"
  empty_cache_key <- checkpointed_genotype_cache_key(
    sample_ids,
    character(),
    design
  )
  base_qc <- checkpointed_preparation_qc(
    raw_dosage_samples,
    overlap_samples,
    sample_ids,
    variant_ids,
    character(),
    plan,
    empty_cache_key
  )
  if (length(sample_ids) == 0L) {
    return(list(
      skip_reason = "TOO_FEW_ALIGNED_SAMPLES",
      genotype = NULL,
      variant_ids = character(),
      design_qr = NULL,
      sample_ids = sample_ids,
      design = design,
      cache_key = empty_cache_key,
      qc = base_qc
    ))
  }
  design_qr <- qr(
    design,
    tol = checkpointed_qr_tolerance,
    LAPACK = FALSE
  )
  if (length(sample_ids) <= design_qr$rank) {
    return(list(
      skip_reason = "TOO_FEW_ALIGNED_SAMPLES",
      genotype = NULL,
      variant_ids = character(),
      design_qr = design_qr,
      sample_ids = sample_ids,
      design = design,
      cache_key = empty_cache_key,
      qc = base_qc
    ))
  }

  genotype[!is.finite(genotype)] <- NA_real_
  imputation_means <- rowMeans(genotype, na.rm = TRUE)
  finite_mean <- is.finite(imputation_means)
  if (!any(finite_mean)) {
    return(list(
      skip_reason = "NO_USABLE_VARIANTS",
      genotype = NULL,
      variant_ids = character(),
      design_qr = design_qr,
      sample_ids = sample_ids,
      design = design,
      cache_key = empty_cache_key,
      qc = base_qc
    ))
  }
  genotype <- genotype[finite_mean, , drop = FALSE]
  retained_variant_ids <- variant_ids[finite_mean]
  imputation_means <- imputation_means[finite_mean]
  missing_cells <- which(is.na(genotype), arr.ind = TRUE)
  if (nrow(missing_cells) > 0L) {
    genotype[missing_cells] <- imputation_means[missing_cells[, "row"]]
  }
  if (!any(genotype > 0)) {
    return(list(
      skip_reason = "NO_ALTERNATIVE_ALLELE",
      genotype = NULL,
      variant_ids = character(),
      design_qr = design_qr,
      sample_ids = sample_ids,
      design = design,
      cache_key = empty_cache_key,
      qc = base_qc
    ))
  }

  variable_variants <- apply(genotype, 1L, function(values) {
    variance <- stats::var(values)
    is.finite(variance) && variance > 0
  })
  if (!any(variable_variants)) {
    return(list(
      skip_reason = "NO_USABLE_VARIANTS",
      genotype = NULL,
      variant_ids = character(),
      design_qr = design_qr,
      sample_ids = sample_ids,
      design = design,
      cache_key = empty_cache_key,
      qc = base_qc
    ))
  }
  genotype <- genotype[variable_variants, , drop = FALSE]
  retained_variant_ids <- retained_variant_ids[variable_variants]

  genotype_samples_by_variants <- t(genotype)
  colnames(genotype_samples_by_variants) <- retained_variant_ids
  genotype_residual <- qr.resid(design_qr, genotype_samples_by_variants)
  residual_variable <- vapply(
    seq_len(ncol(genotype_residual)),
    function(column_index) {
      checkpointed_has_numerical_variation(
        genotype_residual[, column_index],
        genotype_samples_by_variants[, column_index]
      )
    },
    logical(1L)
  )
  if (!any(residual_variable)) {
    return(list(
      skip_reason = "NO_USABLE_VARIANTS",
      genotype = NULL,
      variant_ids = character(),
      design_qr = design_qr,
      sample_ids = sample_ids,
      design = design,
      cache_key = empty_cache_key,
      qc = base_qc
    ))
  }
  genotype_residual <- genotype_residual[, residual_variable, drop = FALSE]
  retained_variant_ids <- retained_variant_ids[residual_variable]
  cache_key <- checkpointed_genotype_cache_key(
    sample_ids,
    retained_variant_ids,
    design
  )
  qc <- checkpointed_preparation_qc(
    raw_dosage_samples,
    overlap_samples,
    sample_ids,
    variant_ids,
    retained_variant_ids,
    plan,
    cache_key
  )
  list(
    skip_reason = NULL,
    genotype = genotype_residual,
    variant_ids = retained_variant_ids,
    design_qr = design_qr,
    sample_ids = sample_ids,
    design = design,
    cache_key = cache_key,
    qc = qc
  )
}

new_checkpointed_window_genotype_cache <- function(
    genotype,
    variant_ids,
    cache_plan,
    raw_dosage_samples,
    overlap_samples
) {
  state <- new.env(parent = emptyenv())
  state$raw_genotype <- genotype
  genotype <- NULL
  state$entries <- list()
  state$prepared_keys <- character()
  state$preparation_counts <- integer()
  plans <- cache_plan$phenotypes
  all_keys <- unique(purrr::map_chr(plans, "preparation_key"))
  last_use <- purrr::map_int(
    all_keys,
    function(key) {
      max(purrr::map_int(
        purrr::keep(plans, ~ identical(.x$preparation_key, key)),
        "processing_index"
      ))
    }
  ) |>
    stats::setNames(all_keys)

  get_entry <- function(phenotype_id) {
    plan <- plans[[phenotype_id]]
    if (is.null(plan)) {
      stop("Cache plan does not contain phenotype: ", phenotype_id, call. = FALSE)
    }
    key <- plan$preparation_key
    if (is.null(state$entries[[key]])) {
      if (is.null(state$raw_genotype)) {
        stop("Raw genotype was released before cache preparation.", call. = FALSE)
      }
      state$entries[[key]] <- prepare_checkpointed_window_genotype(
        state$raw_genotype,
        variant_ids,
        plan,
        raw_dosage_samples,
        overlap_samples
      )
      state$prepared_keys <- c(state$prepared_keys, key)
      current_count <- unname(state$preparation_counts[key])
      if (length(current_count) == 0L || is.na(current_count)) {
        current_count <- 0L
      }
      state$preparation_counts[[key]] <- current_count + 1L
      if (setequal(state$prepared_keys, all_keys)) {
        state$raw_genotype <- NULL
      }
    }
    state$entries[[key]]
  }

  release_entry <- function(phenotype_id) {
    plan <- plans[[phenotype_id]]
    if (is.null(plan)) {
      return(invisible(FALSE))
    }
    key <- plan$preparation_key
    if (identical(plan$processing_index, last_use[[key]])) {
      state$entries[[key]] <- NULL
    }
    invisible(TRUE)
  }

  list(
    get = get_entry,
    release = release_entry,
    metrics = function() {
      list(
        total_preparations = as.integer(sum(state$preparation_counts)),
        preparations_by_key = state$preparation_counts,
        resident_entries = as.integer(length(state$entries)),
        raw_genotype_released = is.null(state$raw_genotype)
      )
    }
  )
}

fit_prepared_checkpointed_window_phenotype <- function(
    prepared,
    phenotype,
    settings
) {
  if (!is.null(prepared$skip_reason)) {
    return(checkpointed_skipped_fit(prepared$skip_reason, prepared$qc))
  }
  phenotype <- if (!is.null(names(phenotype))) {
    as.numeric(phenotype[prepared$sample_ids])
  } else {
    as.numeric(phenotype)
  }
  if (length(phenotype) != length(prepared$sample_ids) ||
      any(!is.finite(phenotype))) {
    stop("Prepared phenotype values must match the cache samples.", call. = FALSE)
  }
  phenotype_variance <- stats::var(phenotype)
  if (!is.finite(phenotype_variance) || phenotype_variance <= 0) {
    return(checkpointed_skipped_fit("ZERO_PHENOTYPE_VARIANCE", prepared$qc))
  }
  transformed_phenotype <- rank_inverse_normal(phenotype)
  phenotype_residual <- qr.resid(prepared$design_qr, transformed_phenotype)
  if (!checkpointed_has_numerical_variation(
    phenotype_residual,
    transformed_phenotype
  )) {
    return(checkpointed_skipped_fit("ZERO_PHENOTYPE_VARIANCE", prepared$qc))
  }

  fit <- susieR::susie(
    prepared$genotype,
    phenotype_residual,
    L = settings$L,
    estimate_residual_variance = settings$estimate_residual_variance,
    estimate_prior_variance = settings$estimate_prior_variance,
    scaled_prior_variance = settings$scaled_prior_variance,
    compute_univariate_zscore = settings$compute_univariate_zscore,
    min_abs_corr = settings$min_abs_corr,
    verbose = TRUE
  )
  fit$variant_id <- prepared$variant_ids
  validate_checkpointed_susie_fit(fit)
  list(
    status = if (isTRUE(fit$converged)) "CONVERGED" else "NONCONVERGED",
    fit = fit,
    reason = NULL,
    qc = prepared$qc
  )
}

fit_checkpointed_window_phenotype <- function(
    genotype,
    phenotype,
    covariates,
    variant_ids,
    settings
) {
  genotype <- as.matrix(genotype)
  storage.mode(genotype) <- "double"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  if (nrow(genotype) != length(variant_ids) ||
      ncol(genotype) != length(phenotype) ||
      nrow(covariates) != length(phenotype)) {
    stop("Genotype, phenotype, and covariate dimensions must align.", call. = FALSE)
  }
  sample_ids <- colnames(genotype) %||% names(phenotype) %||%
    rownames(covariates) %||% sprintf("sample_%08d", seq_len(ncol(genotype)))
  colnames(genotype) <- sample_ids
  if (!is.null(names(phenotype))) {
    phenotype <- phenotype[sample_ids]
  } else {
    names(phenotype) <- sample_ids
  }
  if (!is.null(rownames(covariates)) && all(sample_ids %in% rownames(covariates))) {
    covariates <- covariates[sample_ids, , drop = FALSE]
  }
  complete_samples <- is.finite(phenotype)
  if (ncol(covariates) > 0L) {
    complete_samples <- complete_samples &
      apply(covariates, 1L, function(values) all(is.finite(values)))
  }
  retained_sample_ids <- sample_ids[complete_samples]
  reduced_covariates <- checkpointed_reduce_covariate_matrix(
    covariates[complete_samples, , drop = FALSE]
  )
  design <- cbind(intercept = 1, reduced_covariates)
  rownames(design) <- retained_sample_ids
  plan <- list(
    retained_sample_ids = retained_sample_ids,
    design = design,
    design_id = checkpointed_design_identity(retained_sample_ids, design),
    sample_mask_id = checkpointed_sample_mask_id(retained_sample_ids),
    removed_covariates = attr(reduced_covariates, "removed_covariates") %||%
      character()
  )
  prepared <- prepare_checkpointed_window_genotype(
    genotype,
    variant_ids,
    plan,
    raw_dosage_samples = length(sample_ids),
    overlap_samples = length(sample_ids)
  )
  fit_prepared_checkpointed_window_phenotype(prepared, phenotype, settings)
}

checkpointed_validate_susie_tables <- function(tables) {
  expected <- list(
    credible_sets = names(InitEmptyInCSVariantDf()),
    lbf_variable = names(InitEmptyLbfDf()),
    full_susie = names(InitEmptyVariantDf())
  )
  if (!identical(names(tables), names(expected))) {
    stop("Checkpointed SuSiE tables must use the required names.", call. = FALSE)
  }
  purrr::walk2(tables, expected, function(table, columns) {
    if (!inherits(table, "data.frame") || !identical(names(table), columns)) {
      stop("Checkpointed SuSiE table columns do not match the output schema.", call. = FALSE)
    }
  })
  prototypes <- list(
    credible_sets = InitEmptyInCSVariantDf(),
    lbf_variable = InitEmptyLbfDf(),
    full_susie = InitEmptyVariantDf()
  )
  purrr::walk2(tables, prototypes, function(table, prototype) {
    table_types <- vapply(table, typeof, character(1L))
    prototype_types <- vapply(prototype, typeof, character(1L))
    if (!identical(table_types, prototype_types)) {
      stop(
        "Checkpointed SuSiE table column types do not match the output schema.",
        call. = FALSE
      )
    }
  })
  tables
}

checkpointed_empty_susie_tables <- function() {
  checkpointed_validate_susie_tables(list(
    credible_sets = InitEmptyInCSVariantDf(),
    lbf_variable = InitEmptyLbfDf(),
    full_susie = InitEmptyVariantDf()
  ))
}

checkpointed_call_legacy_helper <- function(helper, ...) {
  helper_environment <- new.env(parent = environment(helper))
  helper_environment$`%>%` <- magrittr::`%>%`
  environment(helper) <- helper_environment
  helper(...)
}

format_checkpointed_susie_tables <- function(result, phenotype_record) {
  if (identical(result$status, "SKIPPED")) {
    return(checkpointed_empty_susie_tables())
  }
  validate_checkpointed_susie_fit(result$fit)
  extraction_fit <- result$fit
  if (is.numeric(extraction_fit$z) && is.null(dim(extraction_fit$z))) {
    extraction_fit$z <- matrix(
      extraction_fit$z,
      ncol = 1L,
      dimnames = list(extraction_fit$variant_id, "z")
    )
  }
  if (is.matrix(extraction_fit$z)) {
    rownames(extraction_fit$z) <- extraction_fit$variant_id
  }
  matrix_fields <- c("alpha", "mu", "mu2", "lbf_variable")
  for (field in matrix_fields) {
    if (!is.null(extraction_fit[[field]]) && is.matrix(extraction_fit[[field]])) {
      colnames(extraction_fit[[field]]) <- extraction_fit$variant_id
    }
  }
  extracted <- checkpointed_call_legacy_helper(extractResults, extraction_fit)
  if (is.null(extracted$variant_df) || nrow(extracted$variant_df) == 0L ||
      is.null(extracted$cs_df) || nrow(extracted$cs_df) == 0L) {
    return(checkpointed_empty_susie_tables())
  }

  phenotype_id <- phenotype_record$phenotype_id[[1L]]
  window_id <- phenotype_record$window_id[[1L]]
  variant_df <- extracted$variant_df |>
    dplyr::mutate(
      phenotype_id = phenotype_id,
      region = window_id
    ) |>
    tidyr::separate(
      .data$variant_id,
      into = c("chr", "pos", "ref", "alt"),
      sep = "_",
      remove = FALSE
    ) |>
    dplyr::mutate(
      chr = stringr::str_remove(.data$chr, "^chr"),
      cs_index = .data$cs_id,
      cs_id = dplyr::if_else(
        is.na(.data$cs_index),
        NA_character_,
        paste(.data$phenotype_id, .data$cs_index, sep = "_")
      )
    )
  pip_df <- checkpointed_call_legacy_helper(extractPipsFromVariantDf, variant_df)
  variant_df <- variant_df |>
    dplyr::left_join(
      pip_df,
      by = c("phenotype_id", "variant_id", "cs_id", "cs_index")
    )

  lbf_variable <- extracted$lbf_df |>
    dplyr::mutate(
      molecular_trait_id = phenotype_id,
      region = window_id
    ) |>
    tidyr::separate(
      .data$variant_id,
      into = c("chromosome", "position", "ref", "alt"),
      sep = "_",
      remove = FALSE
    ) |>
    dplyr::transmute(
      molecular_trait_id = .data$molecular_trait_id,
      region = .data$region,
      variant = .data$variant_id,
      chromosome = stringr::str_remove(.data$chromosome, "^chr"),
      position = as.integer(.data$position),
      dplyr::across(dplyr::all_of(paste0("lbf_variable", seq_len(10L))))
    )

  credible_sets <- variant_df |>
    dplyr::filter(!is.na(.data$cs_index) & !.data$low_purity) |>
    dplyr::transmute(
      molecular_trait_id = .data$phenotype_id,
      variant = .data$variant_id,
      chromosome = .data$chr,
      position = as.integer(.data$pos),
      ref = .data$ref,
      alt = .data$alt,
      cs_id = .data$cs_id,
      cs_index = .data$cs_index,
      region = .data$region,
      pip = .data$pip,
      z = .data$z,
      cs_min_r2 = .data$cs_min_r2,
      cs_avg_r2 = .data$cs_avg_r2,
      cs_size = as.integer(.data$cs_size),
      posterior_mean = .data$posterior_mean,
      posterior_sd = .data$posterior_sd,
      cs_log10bf = .data$cs_log10bf
    )

  component_columns <- c(
    paste0("alpha", seq_len(10L)),
    paste0("mu_", seq_len(10L)),
    paste0("mu2_", seq_len(10L))
  )
  full_susie_prefix <- variant_df |>
    dplyr::transmute(
      molecular_trait_id = .data$phenotype_id,
      variant = .data$variant_id,
      chromosome = .data$chr,
      position = as.integer(.data$pos),
      ref = .data$ref,
      alt = .data$alt,
      cs_id = .data$cs_id,
      cs_index = .data$cs_index,
      low_purity = .data$low_purity,
      region = .data$region,
      pip = .data$pip,
      z = .data$z,
      posterior_mean = .data$posterior_mean,
      posterior_sd = .data$posterior_sd,
      X_column_scale_factors = .data$X_column_scale_factors
    )
  full_susie <- dplyr::bind_cols(
    full_susie_prefix,
    dplyr::select(variant_df, dplyr::all_of(component_columns))
  )
  checkpointed_validate_susie_tables(list(
    credible_sets = credible_sets,
    lbf_variable = lbf_variable,
    full_susie = full_susie
  ))
}

write_checkpointed_susie_tables <- function(tables, output_dir) {
  tables <- checkpointed_validate_susie_tables(tables)
  required <- c("credible_sets", "lbf_variable", "full_susie")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- stats::setNames(
    file.path(output_dir, paste0(required, ".parquet")),
    required
  )
  purrr::walk2(tables, paths, arrow::write_parquet)
  as.list(paths)
}

checkpointed_timestamp <- function(time = Sys.time()) {
  format(
    as.POSIXct(time, tz = "UTC"),
    "%Y-%m-%dT%H:%M:%SZ",
    tz = "UTC"
  )
}

write_local_phenotype_artifacts <- function(result, phenotype_record, directory) {
  if (!is.list(result) || !identical(result$status %in% checkpoint_terminal_statuses, TRUE)) {
    stop("Checkpoint fit result has an invalid status.", call. = FALSE)
  }
  if (!identical(result$status, "SKIPPED")) {
    validate_checkpointed_susie_fit(result$fit)
  }

  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  tables <- format_checkpointed_susie_tables(result, phenotype_record)
  table_paths <- write_checkpointed_susie_tables(tables, directory)
  fit_rds <- NULL
  if (!identical(result$status, "SKIPPED")) {
    fit_rds <- file.path(directory, "susie_fit.rds")
    saveRDS(result$fit, fit_rds)
  }

  c(list(fit_rds = fit_rds), table_paths)
}

checkpointed_payload_record <- function(local_path, relative_path) {
  if (!is.character(local_path) || length(local_path) != 1L ||
      !file.exists(local_path)) {
    stop("Local checkpoint payload is absent.", call. = FALSE)
  }
  list(
    path = relative_path,
    sha256 = sha256_file(local_path),
    bytes = as.numeric(unname(file.info(local_path)$size))
  )
}

checkpointed_package_versions <- function() {
  list(
    R = R.version$version.string,
    arrow = as.character(utils::packageVersion("arrow")),
    susieR = as.character(utils::packageVersion("susieR"))
  )
}

checkpointed_runtime_provenance <- function(config) {
  source_paths <- config$source_paths
  if (is.null(source_paths)) {
    helper_root <- get0(
      "FunctionPath",
      envir = .GlobalEnv,
      inherits = TRUE,
      ifnotfound = file.path(dirname(dirname(config$wrapper_path)), "utils")
    )
    source_paths <- c(
      runner_wrapper = config$wrapper_path,
      checkpointed_functions = file.path(
        helper_root,
        "CheckpointedWindowSusieFunctions.R"
      ),
      checkpoint_store = file.path(helper_root, "CheckpointStore.R")
    )
  }
  required_names <- c(
    "runner_wrapper",
    "checkpointed_functions",
    "checkpoint_store"
  )
  if (!identical(names(source_paths), required_names) ||
      any(!file.exists(source_paths))) {
    stop(
      "Runtime source paths must name the runner wrapper and both checkpoint helpers.",
      call. = FALSE
    )
  }
  runner_image <- config$runtime_image %||%
    Sys.getenv(
      "CHECKPOINTED_SUSIE_RUNTIME_IMAGE",
      unset = "local-unpackaged-runtime"
    )
  base_image_digest <- config$base_image_digest %||%
    Sys.getenv("CHECKPOINTED_SUSIE_BASE_IMAGE_DIGEST")
  if (!nzchar(runner_image) || !nzchar(base_image_digest)) {
    stop("Runtime image identities must be nonempty strings.", call. = FALSE)
  }
  list(
    source_hashes = purrr::map_chr(source_paths, sha256_file),
    runtime = list(
      runner_image = runner_image,
      base_image_digest = base_image_digest
    )
  )
}

checkpointed_complete_fit_qc <- function(
    fit_result,
    plan,
    model_inputs
) {
  qc <- fit_result$qc %||% list()
  retained_variant_ids <- if (identical(fit_result$status, "SKIPPED")) {
    character()
  } else {
    fit_result$fit$variant_id
  }
  defaults <- list(
    input_samples = as.integer(model_inputs$overlap_samples),
    retained_samples = as.integer(length(plan$retained_sample_ids)),
    raw_dosage_samples = as.integer(model_inputs$raw_dosage_samples),
    overlap_samples = as.integer(model_inputs$overlap_samples),
    phenotype_retained_samples = as.integer(length(plan$retained_sample_ids)),
    input_variants = as.integer(nrow(model_inputs$variant_info)),
    retained_variants = as.integer(length(retained_variant_ids)),
    removed_covariates = as.character(plan$removed_covariates),
    removed_variant_ids = model_inputs$variant_info$variant_id[
      !model_inputs$variant_info$variant_id %in% retained_variant_ids
    ],
    design_id = plan$design_id,
    sample_mask_id = plan$sample_mask_id,
    cache_key = checkpointed_genotype_cache_key(
      plan$retained_sample_ids,
      retained_variant_ids,
      plan$design
    )
  )
  for (field in names(defaults)) {
    if (is.null(qc[[field]])) {
      qc[[field]] <- defaults[[field]]
    }
  }
  fit_result$qc <- qc
  fit_result
}

checkpointed_cleanup_phenotype_artifacts <- function(local_artifacts) {
  paths <- unlist(local_artifacts, use.names = FALSE)
  paths <- paths[!is.na(paths) & nzchar(paths)]
  if (length(paths) == 0L) {
    return(invisible(TRUE))
  }
  directories <- unique(dirname(paths))
  if (length(directories) != 1L ||
      !startsWith(basename(directories), "checkpointed-phenotype-")) {
    stop("Phenotype artifact cleanup path is not a managed temporary directory.", call. = FALSE)
  }
  unlink(directories, recursive = TRUE)
  invisible(!dir.exists(directories))
}

build_fit_manifest <- function(
    analysis_id,
    window_id,
    phenotype_record,
    fit_result,
    local_artifacts,
    input_hashes,
    settings,
    provenance
) {
  status <- fit_result$status
  if (!status %in% checkpoint_terminal_statuses) {
    stop("Checkpoint fit result has an invalid status.", call. = FALSE)
  }
  if (!is.list(fit_result$qc)) {
    stop("Checkpoint fit result is missing QC counts.", call. = FALSE)
  }
  container_digest <- provenance$runtime$base_image_digest
  if (!nzchar(container_digest)) {
    stop(
      "CHECKPOINTED_SUSIE_BASE_IMAGE_DIGEST must be a nonempty string.",
      call. = FALSE
    )
  }

  fit_sha256 <- if (identical(status, "SKIPPED")) {
    NULL
  } else {
    sha256_file(local_artifacts$fit_rds)
  }
  paths <- phenotype_checkpoint_paths(
    analysis_id,
    window_id,
    phenotype_record$phenotype_key[[1L]],
    fit_sha256
  )
  artifact_mappings <- checkpoint_artifact_mappings(status)
  payloads <- purrr::map2(
    names(artifact_mappings),
    unname(artifact_mappings),
    function(artifact_name, payload_name) {
      checkpointed_payload_record(
        local_artifacts[[artifact_name]],
        paths[[artifact_name]]
      )
    }
  ) |>
    stats::setNames(unname(artifact_mappings))

  manifest <- list(
    schema_version = checkpoint_schema_version,
    analysis_id = analysis_id,
    window_id = window_id,
    phenotype_id = phenotype_record$phenotype_id[[1L]],
    phenotype_key = phenotype_record$phenotype_key[[1L]],
    modality = phenotype_record$modality[[1L]],
    processing_index = as.integer(phenotype_record$processing_index[[1L]]),
    p_value = as.numeric(phenotype_record$p_value[[1L]]),
    status = status,
    converged = if (identical(status, "SKIPPED")) NULL else fit_result$fit$converged,
    exclusion_reason = if (identical(status, "SKIPPED")) fit_result$reason else NULL,
    input_hashes = as.list(input_hashes),
    settings = settings,
    container_digest = container_digest,
    package_versions = checkpointed_package_versions(),
    source_hashes = as.list(provenance$source_hashes),
    runtime = provenance$runtime,
    counts = fit_result$qc[c(
      "input_samples",
      "retained_samples",
      "raw_dosage_samples",
      "overlap_samples",
      "phenotype_retained_samples",
      "input_variants",
      "retained_variants"
    )],
    qc = fit_result$qc[c(
      "removed_covariates",
      "removed_variant_ids",
      "design_id",
      "sample_mask_id",
      "cache_key"
    )],
    started_at = fit_result$started_at %||% checkpointed_timestamp(),
    completed_at = fit_result$completed_at %||% checkpointed_timestamp(),
    payloads = payloads
  )
  validate_fit_manifest_structure(manifest, phenotype_record)
  manifest
}

checkpointed_download_payload <- function(store, payload, local_path) {
  validate_checkpoint_payload_record(payload, basename(local_path))
  if (!store$object_exists(payload$path)) {
    stop("Checkpoint payload is absent: ", payload$path, call. = FALSE)
  }
  store$download(payload$path, local_path)
  if (!identical(sha256_file(local_path), payload$sha256)) {
    stop("Checkpoint payload checksum does not match: ", payload$path, call. = FALSE)
  }
  if (!identical(
    as.numeric(unname(file.info(local_path)$size)),
    as.numeric(payload$bytes)
  )) {
    stop("Checkpoint payload byte size does not match: ", payload$path, call. = FALSE)
  }
  local_path
}

hydrate_checkpointed_fit_manifests <- function(store, ordered_manifest) {
  purrr::map(
    seq_len(nrow(ordered_manifest)),
    ~ read_valid_boundary_manifest(
      store,
      ordered_manifest[.x, , drop = FALSE]
    )
  )
}

assemble_checkpointed_window_outputs <- function(
    store,
    committed_records,
    output_dir
) {
  if (!is.list(committed_records) || length(committed_records) == 0L) {
    stop("Final assembly requires committed phenotype manifests.", call. = FALSE)
  }
  purrr::walk(committed_records, validate_fit_manifest_structure)
  processing_indexes <- purrr::map_int(
    committed_records,
    ~ checkpoint_scalar_index(.x$processing_index, "Checkpoint processing index")
  )
  record_order <- order(processing_indexes)
  committed_records <- committed_records[record_order]
  processing_indexes <- processing_indexes[record_order]
  if (!identical(processing_indexes, seq_along(committed_records) - 1L)) {
    stop("Committed phenotype manifests are not gap-free and zero-based.", call. = FALSE)
  }

  fit_index <- purrr::map_dfr(committed_records, function(fit_manifest) {
    tibble::tibble(
      window_id = fit_manifest$window_id,
      processing_index = as.integer(fit_manifest$processing_index),
      phenotype_id = fit_manifest$phenotype_id,
      phenotype_key = fit_manifest$phenotype_key,
      modality = fit_manifest$modality,
      p_value = as.numeric(fit_manifest$p_value),
      status = fit_manifest$status,
      converged = fit_manifest$converged %||% NA,
      exclusion_reason = fit_manifest$exclusion_reason %||% NA_character_,
      fit_manifest_uri = fit_manifest$fit_manifest_uri,
      susie_fit_uri = if (identical(fit_manifest$status, "SKIPPED")) {
        NA_character_
      } else {
        fit_manifest$payloads$susie_fit$uri
      },
      susie_fit_sha256 = if (identical(fit_manifest$status, "SKIPPED")) {
        NA_character_
      } else {
        fit_manifest$payloads$susie_fit$sha256
      }
    )
  })

  parquet_names <- c("credible_sets", "lbf_variable", "full_susie")
  combined_tables <- purrr::map(parquet_names, function(payload_name) {
    phenotype_tables <- purrr::map(
      committed_records,
      function(fit_manifest) {
        local_path <- tempfile(
          paste0("assembled-", payload_name, "-"),
          fileext = ".parquet"
        )
        on.exit(unlink(local_path), add = TRUE)
        checkpointed_download_payload(
          store,
          fit_manifest$payloads[[payload_name]],
          local_path
        )
        arrow::read_parquet(local_path, as_data_frame = TRUE)
      }
    )
    dplyr::bind_rows(phenotype_tables)
  }) |>
    stats::setNames(parquet_names)
  combined_tables <- checkpointed_validate_susie_tables(combined_tables)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  fit_index_path <- file.path(output_dir, "window_fit_index.tsv")
  readr::write_tsv(fit_index, fit_index_path, na = "")
  table_paths <- write_checkpointed_susie_tables(combined_tables, output_dir)
  c(list(fit_index = fit_index_path), table_paths)
}

checkpointed_scientific_paths <- function(config) {
  covariate_hash_names <- paste0(
    "covariate_",
    seq_along(config$covariate_files),
    "_",
    config$covariate_modalities
  )
  paths <- c(
    window_dosage = config$window_dosage,
    window_phenotypes = config$window_phenotypes,
    phenotype_data = config$phenotype_data,
    stats::setNames(config$covariate_files, covariate_hash_names)
  )
  if (!is.null(config$keep_samples)) {
    paths <- c(paths, keep_samples = config$keep_samples)
  }
  paths
}

checkpointed_begin_attempt <- function(window_manifest) {
  if (!is.null(window_manifest$failure)) {
    window_manifest$failure_history <- c(
      window_manifest$failure_history %||% list(),
      list(window_manifest$failure)
    )
  }
  window_manifest["failure"] <- list(NULL)
  window_manifest$attempt <- as.integer(window_manifest$attempt %||% 0L) + 1L
  window_manifest$status <- "RUNNING"
  window_manifest
}

run_checkpointed_window <- function(
    config,
    store = NULL,
    fit_function = fit_checkpointed_window_phenotype,
    interrupt_after_commits = Inf
) {
  window_state <- NULL
  window_paths <- NULL
  run_window_body <- function() {
    raw_manifest <- readr::read_tsv(
      config$window_phenotypes,
      col_types = readr::cols(.default = readr::col_character()),
      show_col_types = FALSE
    )
    ordered_manifest <- validate_window_phenotype_manifest(
      raw_manifest,
      config$window_id,
      phenotype_data = config$phenotype_data
    )
    scientific_paths <- checkpointed_scientific_paths(config)
    input_hashes <- purrr::map_chr(scientific_paths, sha256_file)
    settings <- checkpointed_susie_settings()
    provenance <- checkpointed_runtime_provenance(config)
    analysis_id <- build_checkpoint_analysis_id(
      input_hashes = input_hashes,
      ordered_manifest = ordered_manifest,
      settings = settings,
      container_digest = provenance$runtime$base_image_digest,
      wrapper_hash = provenance$source_hashes[["runner_wrapper"]],
      source_hashes = provenance$source_hashes,
      runtime_image = provenance$runtime$runner_image
    )
    ordered_manifest <- ordered_manifest |>
      dplyr::mutate(analysis_id = analysis_id, .before = 1L)

    store <<- store %||% new_checkpoint_store(config$checkpoint_root)
    window_paths <<- window_checkpoint_paths(analysis_id, config$window_id)
    saved_window_manifest <- read_window_run_manifest_if_present(
      store,
      window_paths$window_manifest
    )
    resume <- resolve_resume_boundary(
      store,
      ordered_manifest,
      saved_window_manifest,
      expected_input_hashes = input_hashes,
      expected_settings = settings,
      expected_source_hashes = provenance$source_hashes,
      expected_runtime = provenance$runtime
    )
    if (is.null(resume$window_manifest)) {
      window_state <<- new_window_run_manifest(
        analysis_id,
        config$window_id,
        ordered_manifest,
        input_hashes,
        settings,
        source_hashes = provenance$source_hashes,
        runtime = provenance$runtime
      )
      window_state$last_committed_index <<- resume$last_committed_index
      window_state$recovered_prefix_last_committed_index <<-
        resume$last_committed_index
      window_state$recovery_history <<- resume$recovery_history %||% list()
    } else {
      window_state <<- resume$window_manifest
    }
    window_state <<- checkpointed_begin_attempt(window_state)
    message("[checkpointed-window] Resume at index ", resume$next_index)
    must_persist_resume <- !is.null(saved_window_manifest) ||
      resume$last_committed_index >= 0L
    if (must_persist_resume) {
      upload_window_run_manifest(
        store,
        window_paths$window_manifest,
        window_state
      )
      message(
        "[checkpointed-window] Persisted RUNNING cursor at index ",
        resume$last_committed_index
      )
    }

    commit_count <- 0L
    indexes <- ordered_manifest$processing_index[
      ordered_manifest$processing_index >= resume$next_index
    ]
    model_inputs <- NULL
    cache_plan <- list(phenotypes = list())
    if (length(indexes) > 0L) {
      model_inputs <- load_checkpointed_window_inputs(config, ordered_manifest)
      remaining_manifest <- ordered_manifest |>
        dplyr::filter(.data$processing_index %in% indexes)
      cache_plan <- build_checkpointed_window_cache_plan(
        model_inputs,
        remaining_manifest
      )
    }
    use_prepared_cache <- identical(
      fit_function,
      fit_checkpointed_window_phenotype
    )
    genotype_cache <- NULL
    if (use_prepared_cache && length(indexes) > 0L) {
      genotype_cache <- new_checkpointed_window_genotype_cache(
        model_inputs$genotype,
        model_inputs$variant_info$variant_id,
        cache_plan,
        raw_dosage_samples = model_inputs$raw_dosage_samples,
        overlap_samples = model_inputs$overlap_samples
      )
      model_inputs$genotype <- NULL
    }
    for (index in indexes) {
      phenotype_record <- ordered_manifest |>
        dplyr::filter(.data$processing_index == index)
      phenotype_id <- phenotype_record$phenotype_id[[1L]]
      plan <- cache_plan$phenotypes[[phenotype_id]]
      message(
        "[checkpointed-window] Fit index ",
        index,
        " phenotype ",
        phenotype_id,
        " modality ",
        phenotype_record$modality[[1L]],
        " samples ",
        length(plan$retained_sample_ids),
        " design ",
        plan$design_id
      )
      started_at <- checkpointed_timestamp()
      fit_result <- if (use_prepared_cache) {
        fit_prepared_checkpointed_window_phenotype(
          genotype_cache$get(phenotype_id),
          model_inputs$phenotypes[[phenotype_id]],
          settings
        )
      } else {
        fit_function(
          genotype = model_inputs$genotype,
          phenotype = model_inputs$phenotypes[[phenotype_id]],
          covariates = applicable_covariate_matrix(
            model_inputs$covariates,
            phenotype_record$modality[[1L]],
            model_inputs$sample_ids
          ),
          variant_ids = model_inputs$variant_info$variant_id,
          settings = settings
        )
      }
      fit_result <- checkpointed_complete_fit_qc(
        fit_result,
        plan,
        model_inputs
      )
      fit_result$started_at <- started_at
      fit_result$completed_at <- checkpointed_timestamp()
      local_artifacts <- write_local_phenotype_artifacts(
        fit_result,
        phenotype_record,
        tempfile("checkpointed-phenotype-")
      )
      fit_manifest <- build_fit_manifest(
        analysis_id,
        config$window_id,
        phenotype_record,
        fit_result,
        local_artifacts,
        input_hashes,
        settings,
        provenance
      )
      fit_sha256 <- if (identical(fit_result$status, "SKIPPED")) {
        NULL
      } else {
        fit_manifest$payloads$susie_fit$sha256
      }
      committed <- commit_phenotype_checkpoint(
        store,
        phenotype_checkpoint_paths(
          analysis_id,
          config$window_id,
          phenotype_record$phenotype_key[[1L]],
          fit_sha256
        ),
        local_artifacts,
        fit_manifest
      )
      window_state <<- advance_window_run_manifest(window_state, committed)
      upload_window_run_manifest(
        store,
        window_paths$window_manifest,
        window_state
      )
      checkpointed_cleanup_phenotype_artifacts(local_artifacts)
      if (use_prepared_cache) {
        genotype_cache$release(phenotype_id)
      }
      message(
        "[checkpointed-window] Committed ",
        phenotype_id,
        " at index ",
        index
      )
      commit_count <- commit_count + 1L
      if (commit_count >= interrupt_after_commits) {
        stop(structure(
          list(message = "Synthetic interruption"),
          class = c("checkpoint_test_interrupt", "error", "condition")
        ))
      }
    }

    hydrated_manifests <- hydrate_checkpointed_fit_manifests(
      store,
      ordered_manifest
    )
    window_state$committed <<- purrr::map(
      hydrated_manifests,
      checkpoint_committed_record
    )
    window_state$last_committed_index <<- nrow(ordered_manifest) - 1L
    final_paths <- assemble_checkpointed_window_outputs(
      store,
      hydrated_manifests,
      config$output_dir
    )
    committed_window_outputs <- commit_window_outputs(
      store,
      analysis_id,
      config$window_id,
      final_paths
    )
    window_state <<- complete_window_run_manifest(
      window_state,
      committed_window_outputs
    )
    local_window_manifest <- write_local_window_run_manifest(
      window_state,
      config$output_dir
    )
    upload_window_run_manifest(
      store,
      window_paths$window_manifest,
      window_state
    )
    message("[checkpointed-window] Complete ", config$window_id)
    c(list(window_manifest = local_window_manifest), final_paths)
  }

  tryCatch(
    run_window_body(),
    error = function(condition) {
      if (inherits(condition, "checkpoint_test_interrupt")) {
        stop(condition)
      }
      if (!is.null(window_state) && !is.null(window_paths)) {
        failed_state <- fail_window_run_manifest(window_state, condition)
        try(
          upload_window_run_manifest(
            store,
            window_paths$window_manifest,
            failed_state
          ),
          silent = TRUE
        )
      }
      stop(condition)
    }
  )
}
