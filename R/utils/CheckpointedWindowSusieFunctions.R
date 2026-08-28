# Canonical manifest and identity helpers for one-window univariate SuSiE.

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

checkpoint_phenotype_key <- function(window_id, modality, phenotype_id) {
  digest::digest(
    paste(window_id, modality, phenotype_id, sep = "\n"),
    algo = "sha256",
    serialize = FALSE
  )
}

validate_window_phenotype_manifest <- function(manifest, window_id) {
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
    wrapper_hash
) {
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
    container_digest = container_digest,
    wrapper_hash = wrapper_hash
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
  if (!is.character(variant_id) || length(variant_id) == 0L ||
      anyNA(variant_id) || any(!nzchar(variant_id))) {
    stop("Checkpoint variant IDs must be nonempty strings.", call. = FALSE)
  }
  if (anyDuplicated(variant_id)) {
    stop("Checkpoint variant IDs must be unique.", call. = FALSE)
  }

  pip <- fit$pip
  if (!is.numeric(pip) || length(pip) != length(variant_id) ||
      anyNA(pip) || any(!is.finite(pip)) || any(pip < 0 | pip > 1)) {
    stop(
      "Checkpoint PIP values must be finite and between zero and one.",
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

read_checkpointed_window_dosage <- function(path) {
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

  imputation_means <- rowMeans(genotype, na.rm = TRUE)
  missing_cells <- which(is.na(genotype), arr.ind = TRUE)
  if (nrow(missing_cells) > 0L) {
    replacement_values <- imputation_means[missing_cells[, "row"]]
    finite_replacements <- is.finite(replacement_values)
    genotype[missing_cells[finite_replacements, , drop = FALSE]] <-
      replacement_values[finite_replacements]
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
    c("chromosome", "CHROM", "#CHROM", "chr", "#chr"),
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
  if (any(!is.finite(sample_by_covariate))) {
    stop("Applicable covariate design must be finite.", call. = FALSE)
  }

  variances <- apply(sample_by_covariate, 2L, stats::var)
  nonconstant <- is.finite(variances) & variances > 0
  candidate <- sample_by_covariate[, nonconstant, drop = FALSE]
  retained <- integer()
  if (ncol(candidate) > 0L) {
    for (column_index in seq_len(ncol(candidate))) {
      trial <- candidate[, c(retained, column_index), drop = FALSE]
      trial_rank <- qr(trial, tol = 1e-10, LAPACK = FALSE)$rank
      if (trial_rank > length(retained)) {
        retained <- c(retained, column_index)
      }
    }
  }
  full_rank <- candidate[, retained, drop = FALSE]
  if (ncol(full_rank) == 0L) {
    full_rank <- matrix(
      numeric(length(sample_ids) * 0L),
      nrow = length(sample_ids),
      dimnames = list(sample_ids, character())
    )
  }
  attr(full_rank, "removed_covariates") <- covariate_ids[
    !covariate_ids %in% colnames(full_rank)
  ]
  attr(full_rank, "input_covariates") <- covariate_ids
  full_rank
}

read_checkpointed_sample_allowlist <- function(path) {
  allowlist <- data.table::fread(
    path,
    check.names = FALSE,
    data.table = FALSE
  )
  if (ncol(allowlist) < 1L) {
    stop("Sample allowlist must contain one sample-ID column.", call. = FALSE)
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
  dosage <- read_checkpointed_window_dosage(config$window_dosage)
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
  available_sample_sets <- c(
    purrr::map(phenotype_data$values, names),
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
  finite_phenotypes <- purrr::map(phenotype_data$values, function(values) {
    is.finite(values[aligned_samples])
  }) |>
    purrr::reduce(`&`)
  aligned_samples <- aligned_samples[finite_phenotypes]
  purrr::walk(
    manifest_modalities,
    ~ applicable_covariate_matrix(covariates, .x, aligned_samples)
  )

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
    sample_ids = aligned_samples
  )
}

rank_inverse_normal <- function(values) {
  stats::qnorm(
    (rank(values, ties.method = "average") - 0.5) / length(values)
  )
}

checkpointed_skipped_fit <- function(reason, qc) {
  list(status = "SKIPPED", fit = NULL, reason = reason, qc = qc)
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
  phenotype <- as.numeric(phenotype)
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  variant_ids <- as.character(variant_ids)

  if (nrow(genotype) != length(variant_ids) ||
      ncol(genotype) != length(phenotype) ||
      nrow(covariates) != length(phenotype)) {
    stop("Genotype, phenotype, and covariate dimensions must align.", call. = FALSE)
  }
  if (anyNA(variant_ids) || any(!nzchar(variant_ids)) || anyDuplicated(variant_ids)) {
    stop("Model variant IDs must be unique nonempty strings.", call. = FALSE)
  }

  input_samples <- length(phenotype)
  input_variants <- length(variant_ids)
  complete_samples <- is.finite(phenotype)
  if (ncol(covariates) > 0L) {
    complete_samples <- complete_samples &
      apply(covariates, 1L, function(values) all(is.finite(values)))
  }
  phenotype <- phenotype[complete_samples]
  genotype <- genotype[, complete_samples, drop = FALSE]
  covariates <- covariates[complete_samples, , drop = FALSE]
  retained_samples <- length(phenotype)
  qc <- list(
    input_samples = as.integer(input_samples),
    retained_samples = as.integer(retained_samples),
    input_variants = as.integer(input_variants),
    retained_variants = 0L,
    removed_variant_ids = character()
  )

  design <- cbind(intercept = 1, covariates)
  design_qr <- qr(design)
  if (retained_samples <= design_qr$rank) {
    return(checkpointed_skipped_fit("TOO_FEW_ALIGNED_SAMPLES", qc))
  }
  phenotype_variance <- stats::var(phenotype)
  if (!is.finite(phenotype_variance) || phenotype_variance <= 0) {
    return(checkpointed_skipped_fit("ZERO_PHENOTYPE_VARIANCE", qc))
  }

  genotype[!is.finite(genotype)] <- NA_real_
  imputation_means <- rowMeans(genotype, na.rm = TRUE)
  finite_mean <- is.finite(imputation_means)
  if (!any(finite_mean)) {
    qc$removed_variant_ids <- variant_ids
    return(checkpointed_skipped_fit("NO_USABLE_VARIANTS", qc))
  }
  genotype <- genotype[finite_mean, , drop = FALSE]
  retained_variant_ids <- variant_ids[finite_mean]
  imputation_means <- imputation_means[finite_mean]
  missing_cells <- which(is.na(genotype), arr.ind = TRUE)
  if (nrow(missing_cells) > 0L) {
    genotype[missing_cells] <- imputation_means[missing_cells[, "row"]]
  }
  if (!any(genotype > 0)) {
    qc$removed_variant_ids <- variant_ids
    return(checkpointed_skipped_fit("NO_ALTERNATIVE_ALLELE", qc))
  }

  variable_variants <- apply(genotype, 1L, function(values) {
    variance <- stats::var(values)
    is.finite(variance) && variance > 0
  })
  if (!any(variable_variants)) {
    qc$removed_variant_ids <- variant_ids
    return(checkpointed_skipped_fit("NO_USABLE_VARIANTS", qc))
  }
  genotype <- genotype[variable_variants, , drop = FALSE]
  retained_variant_ids <- retained_variant_ids[variable_variants]

  phenotype_residual <- qr.resid(
    design_qr,
    rank_inverse_normal(phenotype)
  )
  if (!is.finite(stats::var(phenotype_residual)) ||
      stats::var(phenotype_residual) <= 0) {
    return(checkpointed_skipped_fit("ZERO_PHENOTYPE_VARIANCE", qc))
  }
  genotype_samples_by_variants <- t(genotype)
  colnames(genotype_samples_by_variants) <- retained_variant_ids
  genotype_residual <- qr.resid(design_qr, genotype_samples_by_variants)
  residual_variance <- apply(genotype_residual, 2L, stats::var)
  residual_variable <- is.finite(residual_variance) & residual_variance > 0
  if (!any(residual_variable)) {
    qc$removed_variant_ids <- variant_ids
    return(checkpointed_skipped_fit("NO_USABLE_VARIANTS", qc))
  }
  genotype_residual <- genotype_residual[, residual_variable, drop = FALSE]
  retained_variant_ids <- retained_variant_ids[residual_variable]
  qc$retained_variants <- as.integer(length(retained_variant_ids))
  qc$removed_variant_ids <- variant_ids[!variant_ids %in% retained_variant_ids]

  fit <- susieR::susie(
    genotype_residual,
    phenotype_residual,
    L = settings$L,
    estimate_residual_variance = settings$estimate_residual_variance,
    estimate_prior_variance = settings$estimate_prior_variance,
    scaled_prior_variance = settings$scaled_prior_variance,
    compute_univariate_zscore = settings$compute_univariate_zscore,
    min_abs_corr = settings$min_abs_corr,
    verbose = TRUE
  )
  fit$variant_id <- retained_variant_ids
  validate_checkpointed_susie_fit(fit)
  list(
    status = if (isTRUE(fit$converged)) "CONVERGED" else "NONCONVERGED",
    fit = fit,
    reason = NULL,
    qc = qc
  )
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
