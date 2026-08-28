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
