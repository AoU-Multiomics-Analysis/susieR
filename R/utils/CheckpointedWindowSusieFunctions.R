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
