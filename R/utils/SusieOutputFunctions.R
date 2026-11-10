# Shared output formatting for the legacy and PGEN entry points.
# Keep the native fitted values and append the recorded phenotype scale.
# Optional extractions must come from these same fitted outcomes.
formatSusieOutputs <- function(results, region_df, extracted_results = NULL) {
  outcome_ids <- names(results)
  if (!is.list(results) || (length(results) &&
      (is.null(outcome_ids) || anyNA(outcome_ids) || any(!nzchar(outcome_ids)) ||
       anyDuplicated(outcome_ids)))) {
    stop("results must be a named list with unique nonempty outcome IDs", call. = FALSE)
  }
  if (!all(c("phenotype_id", "region") %in% names(region_df)) ||
      anyNA(region_df$phenotype_id) || anyDuplicated(region_df$phenotype_id) ||
      !all(outcome_ids %in% region_df$phenotype_id)) {
    stop("region_df must map each unique outcome ID to one region", call. = FALSE)
  }
  if (!is.null(extracted_results)) {
    extracted_ids <- names(extracted_results)
    if (!is.list(extracted_results) || (length(extracted_results) &&
        (is.null(extracted_ids) || anyNA(extracted_ids) || any(!nzchar(extracted_ids)) ||
         anyDuplicated(extracted_ids)))) {
      stop("extracted_results must be a named list with unique nonempty outcome IDs", call. = FALSE)
    }
    if (length(extracted_results) != length(results) || !setequal(extracted_ids, outcome_ids)) {
      stop("extracted_results must contain exactly one extraction for each outcome", call. = FALSE)
    }
    required_tables <- c("cs_df", "variant_df", "lbf_df")
    complete <- vapply(extracted_results, function(extracted) {
      is.list(extracted) && all(required_tables %in% names(extracted)) &&
        !anyDuplicated(names(extracted))
    }, logical(1))
    if (!all(complete)) {
      stop("Each extraction must contain cs_df, variant_df, and lbf_df", call. = FALSE)
    }
    valid_tables <- vapply(extracted_results, function(extracted) {
      all(vapply(extracted[required_tables], function(table) {
        is.null(table) || is.data.frame(table)
      }, logical(1)))
    }, logical(1))
    if (!all(valid_tables)) {
      stop("Extracted tables must be data frames or NULL", call. = FALSE)
    }
    # Use fitted-outcome order even when supplied extractions have another order.
    extracted_results <- extracted_results[outcome_ids]
  }
  empty_variant_df <- InitEmptyVariantDf(include_phenotype_residual_sd = TRUE)
  empty_lbf_df <- InitEmptyLbfDf()
  empty_in_cs_variant_df <- InitEmptyInCSVariantDf(
    include_covariate_diagnostics = TRUE, include_phenotype_residual_sd = TRUE
  )
  empty_cs_df <- InitEmptyCS()
  covariate_explained_variants <- dplyr::bind_rows(
    dplyr::tibble(molecular_trait_id = character(), variant = character(),
                  input_scale = double(), residual_sd = double(),
                  residual_sd_ratio = double()),
    purrr::imap_dfr(results, function(fitted, phenotype_id) {
      diagnostic <- fitted$covariate_explained_variants
      if (is.null(diagnostic) || nrow(diagnostic) == 0L) return(NULL)
      diagnostic %>%
        dplyr::transmute(molecular_trait_id = phenotype_id, variant = variant_id,
                         input_scale, residual_sd, residual_sd_ratio)
    })
  )
  if (is.null(extracted_results)) extracted_results <- purrr::map(results, extractResults)
  res <- purrr::transpose(extracted_results)
  variant_df <- purrr::map_df(res$variant_df, identity, .id = "phenotype_id")
  if (nrow(variant_df) > 0L) {
    variant_df <- variant_df %>%
      dplyr::left_join(region_df, by = "phenotype_id") %>%
      tidyr::separate(variant_id, c("chr", "pos", "ref", "alt"), sep = "_", remove = FALSE) %>%
      dplyr::mutate(chr = stringr::str_remove_all(chr, "chr")) %>%
      dplyr::mutate(cs_index = cs_id) %>%
      dplyr::mutate(cs_id = paste(phenotype_id, cs_index, sep = "_"))
  }

  lbf_df_res <- purrr::map_df(res$lbf_df, identity, .id = "phenotype_id")
  if (nrow(lbf_df_res) > 0L) {
    lbf_df <- lbf_df_res %>%
      dplyr::left_join(region_df, by = "phenotype_id") %>%
      tidyr::separate(variant_id, c("chromosome", "position", "ref", "alt"),
                      sep = "_", remove = FALSE) %>%
      dplyr::mutate(chromosome = stringr::str_remove_all(chromosome, "chr")) %>%
      dplyr::select(-ref, -alt) %>%
      dplyr::rename(molecular_trait_id = phenotype_id, variant = variant_id) %>%
      dplyr::select(molecular_trait_id, region, variant, chromosome, position,
                    lbf_variable1:lbf_variable10) %>%
      dplyr::mutate(position = as.integer(position))
  } else {
    lbf_df <- empty_lbf_df
  }

  cs_df <- purrr::map_df(res$cs_df, identity, .id = "phenotype_id")
  if (nrow(cs_df) > 0L) {
    cs_df <- dplyr::left_join(cs_df, region_df, by = "phenotype_id") %>%
      dplyr::mutate(cs_index = cs_id) %>%
      dplyr::mutate(cs_id = paste(phenotype_id, cs_index, sep = "_")) %>%
      dplyr::transmute(molecular_trait_id = phenotype_id, cs_id, cs_index, region,
                       cs_log10bf, cs_avg_r2, cs_min_r2, cs_size, low_purity,
                       has_covariate_explained_variant)
    in_cs_variant_df <- dplyr::filter(variant_df, !is.na(cs_index) & !low_purity) %>%
      dplyr::transmute(molecular_trait_id = phenotype_id, variant = variant_id,
                       chromosome = chr, position = pos, ref, alt, cs_id, cs_index,
                       region, pip, z, cs_min_r2, cs_avg_r2, cs_size, posterior_mean,
                       posterior_sd, cs_log10bf, has_covariate_explained_variant,
                       phenotype_residual_SD)
  } else {
    in_cs_variant_df <- empty_in_cs_variant_df
    cs_df <- empty_cs_df
  }

  if (nrow(variant_df) > 0L) {
    variant_df_transmute <- dplyr::transmute(variant_df,
      molecular_trait_id = phenotype_id, variant = variant_id, chromosome = chr,
      position = pos, ref, alt, cs_id, cs_index, low_purity, region, pip, z,
      posterior_mean, posterior_sd, X_column_scale_factors)
    variant_df <- dplyr::bind_cols(variant_df_transmute,
                                   dplyr::select(variant_df, alpha1:mu2_10,
                                                 phenotype_residual_SD))
  } else {
    variant_df <- empty_variant_df
  }

  # The old no-CS path writes the schema constructors without a process exit.
  if (nrow(variant_df) == 0L && nrow(cs_df) == 0L && nrow(in_cs_variant_df) == 0L) {
    lbf_df <- empty_lbf_df
    variant_df <- empty_variant_df
  } else {
    in_cs_variant_df <- dplyr::mutate(in_cs_variant_df, position = as.integer(position))
  }
  list(in_cs_variant_df = in_cs_variant_df, lbf_df = lbf_df,
       variant_df = variant_df, cs_df = cs_df,
       covariate_explained_variants = covariate_explained_variants)
}

# Write the same five files as susie.R. Return paths for callers to collect.
writeSusieOutputs <- function(results, tables, output_prefix) {
  files <- c(rds = paste0(output_prefix, "_susie.rds"),
             summary = paste0(output_prefix, ".parquet"),
             lbf = paste0(output_prefix, ".lbf_variable.parquet"),
             full = paste0(output_prefix, ".full_susie.parquet"),
             diagnostics = paste0(output_prefix, ".covariate_explained_variants.tsv"))
  saveRDS(results, files[["rds"]])
  readr::write_tsv(tables$covariate_explained_variants, files[["diagnostics"]],
                   quote = "needed")
  arrow::write_parquet(tables$in_cs_variant_df, files[["summary"]])
  arrow::write_parquet(tables$lbf_df, files[["lbf"]])
  arrow::write_parquet(tables$variant_df, files[["full"]])
  files
}
