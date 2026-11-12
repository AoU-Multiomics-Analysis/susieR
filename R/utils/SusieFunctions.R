# Core SuSiE fine-mapping helpers used by susie.R and the CV workflow.

# Filter variants by per-ancestry MAF or by an explicit variant allowlist.
filterMAF <- function(genotype_matrix,
                      ancestry_df,
                      MAF_threshold = 0,
                      variant_list = NULL) {
    # If no variant list is specified, compute MAFs from the dosage matrix per
    # ancestry and keep variants passing in at least one non-admixed ancestry.
    if (is.null(variant_list) && MAF_threshold > 0 && !is.null(ancestry_df))  {
    message('Filtering genotype matrix using MAF from dataset per ancestry')
    MAF_calculations_table <- genotype_matrix %>% 
        t() %>% 
        data.frame() %>% 
        mutate(across(everything(),~case_when(. == -1 ~ NA,TRUE ~ .))) %>% 
        tibble::rownames_to_column('research_id') %>%
        mutate(research_id = as.numeric(research_id)) %>% 
        left_join(ancestry_df,by = 'research_id') %>%
        tibble::column_to_rownames('research_id') %>% 
        filter(ancestry_pred_other != 'oth') %>% 
        group_by(ancestry_pred_other) %>% 
          summarize(
            across(
              everything(),
              ~ sum(., na.rm = TRUE) / (sum(!is.na(.)) * 2)
            ),
            .groups = "drop"
          ) 
     MAF_calculations <- MAF_calculations_table  %>% 
        pivot_longer(!ancestry_pred_other)   %>% 
        mutate(value = case_when(value > .5 ~ 1 - value,TRUE ~ value)) %>% 
        filter(value > MAF_threshold) %>% 
        distinct(name) %>% 
        pull(name)
    filtered_genotype_matrix <- genotype_matrix[MAF_calculations,]
    }else if (!is.null(variant_list)){
        message('Filtering genotype matrix by external  variant set')
        # Variant lists take precedence over MAF filtering.
        variant_df <- ImportVariantList(variant_list)
        genotype_df_variant_names <- data.frame(variants=rownames(genotype_matrix) ) %>% 
                                            mutate(variants = str_replace(variants,'chrchr','chr')) %>% 
                                            pull(variants)
        rownames(genotype_matrix) <- genotype_df_variant_names
        valid_variants <- intersect(rownames(genotype_matrix),variant_df)
        filtered_genotype_matrix <- genotype_matrix[valid_variants,]
    } else {
        message('No filtering parameters supplied, returning original genotype matrix')
        filtered_genotype_matrix <- genotype_matrix
    }

    return(filtered_genotype_matrix) 
}


# Impute missing values (sentinel -1 or NA) with row means.
impute_matrix_rowmean <- function(mat, missing_value = -1) {
  # mat: numeric matrix (variants x samples)
  # Replace sentinel with NA
  mat[mat == missing_value] <- NA_real_
  rmv <- rowMeans(mat, na.rm = TRUE)
  # If a row is all NA, set to 0
  rmv[is.na(rmv)] <- 0
  # Fill NAs by broadcasting row means
  inds_na <- which(is.na(mat), arr.ind = TRUE)
  if (nrow(inds_na) > 0) {
    mat[inds_na] <- rmv[inds_na[,1]]
  }
  mat
}

# Standardize variants by subtracting row means without large transposes.
standardize_rows <- function(mat) {
  rm <- rowMeans(mat)
  # subtract row means (in-place-like by using sweep)
  mat <- mat - rm
  # return (variants x samples) standardized matrix
  mat
}

# Record frequency from complete, uncentered diploid ALT dosages. Generic
# additional features outside this encoding have unknown allele frequency.
cohortMafFromDosages <- function(mat) {
  if (!nrow(mat)) return(numeric())
  frequency <- rowMeans(mat) / 2
  ranges <- matrixStats::rowRanges(mat)
  maf <- pmin(frequency, 1 - frequency)
  known <- is.finite(ranges[, 1L]) & is.finite(ranges[, 2L]) &
    ranges[, 1L] >= 0 & ranges[, 2L] <= 2 & is.finite(maf)
  maf[!known] <- NA_real_
  unname(as.double(maf))
}

normalize_chr <- function(chr) {
  stringr::str_remove(as.character(chr), "^(?i:chr)+")
}

parse_variant_ids <- function(variant_ids) {
  parts <- stringr::str_split_fixed(variant_ids, "_", 4)
  data.frame(
    variant_id = variant_ids,
    chromosome = normalize_chr(parts[, 1]),
    position = suppressWarnings(as.integer(parts[, 2])),
    stringsAsFactors = FALSE
  )
}

writeVariantPositionSummary <- function(summary_df, output_file) {
  if (is.null(output_file)) {
    return(invisible(NULL))
  }

  write_header <- !file.exists(output_file)
  quoted_columns <- which(vapply(summary_df, function(column) {
    is.character(column) && any(grepl('["\t\r\n]', column), na.rm = TRUE)
  }, logical(1)))
  utils::write.table(
    summary_df,
    file = output_file,
    sep = "\t",
    quote = if (length(quoted_columns)) quoted_columns else FALSE,
    qmethod = "double",
    row.names = FALSE,
    col.names = write_header,
    append = !write_header
  )
  invisible(NULL)
}

reportVariantDirectionCounts <- function(variant_ids,
                                         feature_meta,
                                         context = "genotype matrix",
                                         output_file = NULL,
                                         requested_window_start = NA_real_,
                                         requested_window_end = NA_real_) {
  variant_info <- parse_variant_ids(variant_ids)
  feature_meta <- feature_meta[seq_len(min(nrow(feature_meta), 1)), , drop = FALSE]

  if (nrow(feature_meta) == 0) {
    message("Variant position summary skipped: no feature metadata available")
    return(invisible(NULL))
  }

  feature_chr <- normalize_chr(feature_meta$chromosome[[1]])
  feature_pos <- suppressWarnings(as.numeric(feature_meta$phenotype_pos[[1]]))
  feature_id <- feature_meta$phenotype_id[[1]]

  if (is.na(feature_pos)) {
    message(
      "Variant position summary skipped for ",
      feature_id,
      ": feature position is missing or non-numeric"
    )
    return(invisible(NULL))
  }

  same_chr <- variant_info$chromosome == feature_chr & !is.na(variant_info$position)
  signed_distances <- variant_info$position[same_chr] - feature_pos
  abs_distances <- abs(signed_distances)
  upstream_count <- sum(same_chr & variant_info$position < feature_pos)
  downstream_count <- sum(same_chr & variant_info$position > feature_pos)
  at_feature_count <- sum(same_chr & variant_info$position == feature_pos)
  uncounted_count <- sum(!same_chr | is.na(variant_info$position))
  min_signed_distance <- if (length(signed_distances) > 0) min(signed_distances) else NA_real_
  max_signed_distance <- if (length(signed_distances) > 0) max(signed_distances) else NA_real_
  min_abs_distance <- if (length(abs_distances) > 0) min(abs_distances) else NA_real_
  max_abs_distance <- if (length(abs_distances) > 0) max(abs_distances) else NA_real_
  requested_window_size <- if (
    !is.na(requested_window_start) &&
      !is.na(requested_window_end)
  ) {
    requested_window_end - requested_window_start + 1
  } else {
    NA_real_
  }
  summary_df <- data.frame(
    phenotype_id = feature_id,
    context = context,
    feature_chromosome = feature_chr,
    feature_position = feature_pos,
    requested_window_start = requested_window_start,
    requested_window_end = requested_window_end,
    requested_window_size_bp = requested_window_size,
    total_variants_tested = length(variant_ids),
    same_chromosome_parseable_variants = sum(same_chr),
    upstream_variants = upstream_count,
    downstream_variants = downstream_count,
    variants_at_feature_coordinate = at_feature_count,
    uncounted_variants = uncounted_count,
    min_abs_distance_bp = min_abs_distance,
    max_abs_distance_bp = max_abs_distance,
    min_signed_distance_bp = min_signed_distance,
    max_signed_distance_bp = max_signed_distance,
    stringsAsFactors = FALSE
  )

  message(
    "Variant position summary for ",
    feature_id,
    " at ",
    feature_chr,
    ":",
    feature_pos,
    " after loading ",
    context,
    ": ",
    length(variant_ids),
    " total variants tested, ",
    upstream_count,
    " upstream, ",
    downstream_count,
    " downstream, ",
    at_feature_count,
    " at the feature coordinate, ",
    uncounted_count,
    " uncounted due to chromosome mismatch or unparsable position, ",
    "min abs distance ",
    min_abs_distance,
    " bp, max abs distance ",
    max_abs_distance,
    " bp, requested window ",
    requested_window_start,
    "-",
    requested_window_end,
    " (",
    requested_window_size,
    " bp)"
  )
  writeVariantPositionSummary(summary_df, output_file)
  invisible(NULL)
}

merge_phenotype_windows <- function(phenotype_meta, phenotype_ids, cis_distance) {
  windows <- phenotype_meta %>%
    dplyr::filter(phenotype_id %in% phenotype_ids) %>%
    dplyr::transmute(
      chromosome = as.character(chromosome),
      chromosome_norm = normalize_chr(chromosome),
      start = pmax(1, as.numeric(phenotype_pos) - cis_distance),
      end = as.numeric(phenotype_pos) + cis_distance
    ) %>%
    dplyr::arrange(chromosome_norm, start, end)

  if (nrow(windows) == 0) {
    return(windows)
  }

  merged <- list()
  current <- windows[1, , drop = FALSE]
  for (i in seq_len(nrow(windows))[-1]) {
    row <- windows[i, , drop = FALSE]
    if (row$chromosome_norm == current$chromosome_norm && row$start <= current$end) {
      current$end <- max(current$end, row$end)
    } else {
      merged[[length(merged) + 1]] <- current
      current <- row
    }
  }
  merged[[length(merged) + 1]] <- current
  dplyr::bind_rows(merged)
}

# Store the small covariate factor instead of a sample-by-sample matrix.
emptyCovariateExplainedVariants <- function() {
  tibble::tibble(variant_id = character(), input_scale = double(),
                 residual_sd = double(), residual_sd_ratio = double())
}

createCovariateAdjustment <- function(covariates_matrix) {
  state <- new.env(parent = emptyenv())
  state$hat <- NULL
  state$legacy_cache <- NULL
  state$covariate_explained_variants <- emptyCovariateExplainedVariants()
  list(
    covariates = covariates_matrix,
    projection = covariates_matrix %*% solve(crossprod(covariates_matrix)),
    state = state
  )
}

# Inputs must already follow the covariate sample order, as in the dense path.
residualizeWithCovariates <- function(adjustment, values, force_dense = FALSE,
                                      variant_ids = NULL) {
  dense <- force_dense || !is.null(adjustment$state$hat)
  if (dense) {
    if (is.null(adjustment$state$hat)) {
      adjustment$state$hat <- diag(nrow(adjustment$covariates)) -
        adjustment$projection %*% t(adjustment$covariates)
    }
    residuals <- adjustment$state$hat %*% values
    if (is.null(variant_ids)) return(residuals)
  } else {
    residuals <- values - adjustment$projection %*% crossprod(adjustment$covariates, values)
  }
  # SuSiE can amplify rounding errors when it scales an almost fully
  # explained column. Compare column scales without a full squared matrix.
  input_matrix <- if (is.null(dim(values))) matrix(values, ncol = 1L) else values
  input_scale <- pmax(matrixStats::colSds(input_matrix), abs(colMeans(input_matrix)))
  residual_sd <- matrixStats::colSds(residuals)
  if (!dense && any(input_scale > 0 & residual_sd <= 1e-8 * input_scale, na.rm = TRUE)) {
    # Keep the original multiplication dimensions and operation order for the
    # complete input. A column subset can use different BLAS accumulation.
    adjustment$state$hat <- diag(nrow(adjustment$covariates)) -
      adjustment$projection %*% t(adjustment$covariates)
    residuals <- adjustment$state$hat %*% values
    if (!is.null(variant_ids)) residual_sd <- matrixStats::colSds(residuals)
  }
  if (!is.null(variant_ids)) {
    if (length(variant_ids) != ncol(residuals)) {
      stop("Variant IDs must match the genotype columns")
    }
    flagged <- which(input_scale > 0 & residual_sd <= 1e-8 * input_scale)
    adjustment$state$covariate_explained_variants <- tibble::tibble(
      variant_id = as.character(variant_ids[flagged]),
      input_scale = unname(input_scale[flagged]),
      residual_sd = unname(residual_sd[flagged]),
      residual_sd_ratio = unname(residual_sd[flagged] / input_scale[flagged])
    )
  }
  residuals
}

prepareReusableGenotype <- function(phenotype_ids,
                                    se,
                                    genotype_file,
                                    covariates,
                                    cis_distance,
                                    ancestry_data = NULL,
                                    MAF = 0,
                                    variant_list = NULL,
                                    additional_genotypes = NULL,
                                    force_dense = FALSE) {
  phenotype_ids <- unname(unlist(phenotype_ids))
  if (!is.character(genotype_file)) {
    message("Reusable genotype matrix requested but genotype input is already pre-loaded; using per-phenotype path")
    return(NULL)
  }
  if (!is.null(additional_genotypes)) {
    message("Reusable genotype matrix requested but additional genotypes are present; using per-phenotype path")
    return(NULL)
  }

  phenotype_meta <- SummarizedExperiment::rowData(se) %>% as.data.frame()
  merged_windows <- merge_phenotype_windows(phenotype_meta, phenotype_ids, cis_distance)
  if (nrow(merged_windows) != 1) {
    message(
      "Reusable genotype matrix requested, but selected phenotype windows merged into ",
      nrow(merged_windows),
      " regions; using per-phenotype path"
    )
    return(NULL)
  }

  message(
    "Preparing reusable genotype matrix for ",
    merged_windows$chromosome,
    ":",
    merged_windows$start,
    "-",
    merged_windows$end
  )
  first_phenotype <- phenotype_ids[[1]]
  gene_vector <- eQTLUtils::extractPhentypeFromSE(first_phenotype, se, "counts")
  sample_ids <- gene_vector$genotype_id
  covariates_matrix <- cbind(covariates[sample_ids, ], 1)
  adjustment <- createCovariateAdjustment(covariates_matrix)

  genotype_matrix <- eQTLUtils::extractGenotypeMatrixFromDosage(
    chr = merged_windows$chromosome,
    start = merged_windows$start,
    end = merged_windows$end,
    dosage_file = genotype_file
  ) %>%
    filterMAF(ancestry_data, MAF_threshold = MAF, variant_list = variant_list)

  message("Subsetting reusable genotype matrix to individuals in the molecular data")
  gt_matrix <- genotype_matrix[, sample_ids]
  rm(genotype_matrix)
  gc()

  nonzero_idx <- which(rowSums(round(gt_matrix, 0), na.rm = TRUE) != 0)
  if (length(nonzero_idx) == 0) {
    stop("No variant with alt alleles")
  }
  gt_matrix <- gt_matrix[nonzero_idx, , drop = FALSE]
  all_missing <- matrixStats::rowCounts(gt_matrix, value = NA_real_) +
    matrixStats::rowCounts(gt_matrix, value = -1, na.rm = TRUE) == ncol(gt_matrix)
  gt_matrix <- impute_matrix_rowmean(gt_matrix, missing_value = -1)
  variant_names <- rownames(gt_matrix)
  cohort_maf <- cohortMafFromDosages(gt_matrix)
  cohort_maf[all_missing] <- NA_real_

  message("Standardizing reusable genotype matrix")
  gt_std <- standardize_rows(gt_matrix)
  rm(gt_matrix)
  gc()

  message("Computing reusable gt_hat")
  gt_hat <- residualizeWithCovariates(adjustment, t(gt_std), force_dense = force_dense,
                                     variant_ids = variant_names)
  rm(gt_std)
  gc()

  list(
    gt_hat = gt_hat,
    variant_names = variant_names,
    cohort_maf = cohort_maf,
    variant_info = parse_variant_ids(variant_names),
    sample_ids = sample_ids,
    covariate_adjustment = adjustment,
    dense = !is.null(adjustment$state$hat),
    covariate_explained_variants = adjustment$state$covariate_explained_variants,
    phenotype_ids = phenotype_ids
  )
}

#MergeAdditionalGenotypes <- function(GenotypeMatrix,AdditionalGenotypes) {

#GenotypeMatrixSampleList <- as.character(colnames(GenotypeMatrix))
#AdditionalGenotypesSampleList <- as.character(colnames(AdditionalGenotypes))

#GenotypeMatrixSubset <- GenotypeMatrix %>% 
                            #data.frame() %>%
                            #dplyr::rename_with(~str_remove(.,'^X')) %>% 
                            #select(any_of(AdditionalGenotypesSampleList))
#AdditionalGenotypesSubset <- AdditionalGenotypes %>% 
                            #data.frame() %>% 
                            #dplyr::rename_with(~str_remove(.,'^X')) %>% 
                            #select(any_of(GenotypeMatrixSampleList))

#MergedGenotypes <- bind_rows(GenotypeMatrixSubset,
                             #AdditionalGenotypesSubset
                            #) %>% 
                    #data.matrix()
#MergedGenotypes
#}

# Merge ordinary dosage data with optional additional genotype features while
# preserving only samples shared by both matrices.
MergeAdditionalGenotypes <- function(GenotypeMatrix, AdditionalGenotypes) {
  common_samples <- intersect(
    str_remove(colnames(GenotypeMatrix), "^X"),
    str_remove(colnames(AdditionalGenotypes), "^X")
  )

  if (length(common_samples) == 0) {
    stop("No overlapping samples between GenotypeMatrix and AdditionalGenotypes")
  }

  gt_var_ids  <- rownames(GenotypeMatrix)
  add_var_ids <- rownames(AdditionalGenotypes)

  GenotypeMatrixSubset <- GenotypeMatrix %>%
    data.frame(check.names = FALSE) %>%
    dplyr::rename_with(~str_remove(., "^X")) %>%
    dplyr::select(all_of(common_samples))

  AdditionalGenotypesSubset <- AdditionalGenotypes %>%
    data.frame(check.names = FALSE) %>%
    dplyr::rename_with(~str_remove(., "^X")) %>%
    dplyr::select(all_of(common_samples))

  rownames(GenotypeMatrixSubset) <- gt_var_ids
  rownames(AdditionalGenotypesSubset) <- add_var_ids

  MergedGenotypes <- dplyr::bind_rows(
    GenotypeMatrixSubset,
    AdditionalGenotypesSubset
  ) %>%
    data.matrix()

  return(MergedGenotypes)
}

# Keep the original model settings and calculate marginal scores after fitting.
fitSusieModel <- function(X, y, legacy = FALSE,
                          covariate_explained_variants = emptyCovariateExplainedVariants()) {
  message("Fitting SuSiE model")
  postprocess <- if (legacy) NULL else function(fitted, stats) {
    message("Calculating marginal z-scores")
    fitted$z <- computeMarginalZScores(X, y, .x_stats = stats)
    fitted
  }
  fitted <- runSusieWithFitOptimizations(X, y,
                        L = 10,
                        estimate_residual_variance = TRUE,
                        estimate_prior_variance = TRUE,
                        scaled_prior_variance = 0.1,
                        verbose = TRUE,
                        compute_univariate_zscore = legacy,
                        min_abs_corr = 0.5,
                        .postprocess = postprocess)
  fitted$covariate_explained_variants <- covariate_explained_variants
  fitted$phenotype_residual_SD <- stats::sd(drop(y))
  fitted
}

# Fine-map one phenotype by loading local cis dosages, residualizing expression
# and genotypes against covariates, then fitting susieR.
finemapPhenotype <- function(phenotype_id, se, genotype_file, covariates, cis_distance,ancestry_data = NULL,MAF = 0,variant_list = NULL,additional_genotypes = NULL,reusable_genotype = NULL,variant_report_file = NULL,timing_callback = NULL,genotype_imputer = NULL,genotypes_imputed = FALSE){
  if (!is.null(timing_callback) && !is.function(timing_callback)) {
    stop("timing_callback must be a function or NULL", call. = FALSE)
  }
  timing_event <- function(stage, event) {
    if (!is.null(timing_callback)) timing_callback(stage, event)
    invisible(NULL)
  }
  timing_stage <- "data_preparation"
  timing_event(timing_stage, "start")
  on.exit({
    if (!is.null(timing_stage)) timing_event(timing_stage, "failed")
  }, add = TRUE)
  timed_fit <- function(...) {
    timing_event("data_preparation", "complete")
    timing_stage <<- "fine_mapping"
    timing_event(timing_stage, "start")
    fitted <- fitSusieModel(...)
    timing_event(timing_stage, "complete")
    timing_stage <<- NULL
    fitted
  }
  message("Processing phenotype: ", phenotype_id)
  
  # Extract and rank-normalize the target phenotype from the SummarizedExperiment.
  gene_vector = eQTLUtils::extractPhentypeFromSE(phenotype_id, se, "counts") %>%
    dplyr::mutate(phenotype_value_std = qnorm((rank(phenotype_value, na.last = "keep") - 0.5) / sum(!is.na(phenotype_value))))
  selected_phenotype = phenotype_id
  gene_meta = dplyr::filter(SummarizedExperiment::rowData(se) %>% as.data.frame(), phenotype_id == selected_phenotype)

  if (!is.null(reusable_genotype)) {
    message("Using reusable residualized genotype matrix")
    cache_state <- reusable_genotype$covariate_adjustment$state
    if (!is.null(cache_state$legacy_cache)) {
      reusable_genotype <- cache_state$legacy_cache
    }
    sample_ids <- reusable_genotype$sample_ids
    gene_vector <- gene_vector %>%
      dplyr::filter(genotype_id %in% sample_ids)
    sample_match <- match(sample_ids, gene_vector$genotype_id)
    if (any(is.na(sample_match))) {
      stop("Reusable genotype matrix sample set does not match phenotype samples")
    }
    gene_vector <- gene_vector[sample_match, , drop = FALSE]
    expression_vector <- residualizeWithCovariates(reusable_genotype$covariate_adjustment,
                                                  gene_vector$phenotype_value_std)
    if (!reusable_genotype$dense &&
        !is.null(reusable_genotype$covariate_adjustment$state$hat)) {
      # A later response can require legacy arithmetic after a fast genotype
      # cache was made. Rebuild the same full merged window, then subset it.
      message("Rebuilding the reusable genotype matrix with the legacy covariate adjustment")
      legacy_cache <- prepareReusableGenotype(
        reusable_genotype$phenotype_ids, se, genotype_file, covariates,
        cis_distance, ancestry_data = ancestry_data, MAF = MAF,
        variant_list = variant_list, additional_genotypes = additional_genotypes,
        force_dense = TRUE
      )
      if (!identical(legacy_cache$sample_ids, reusable_genotype$sample_ids) ||
          !identical(legacy_cache$variant_names, reusable_genotype$variant_names)) {
        stop("Legacy reusable genotype rebuild changed the sample or variant order")
      }
      cache_state$legacy_cache <- legacy_cache
      cache_state$hat <- NULL
      reusable_genotype <- legacy_cache
      expression_vector <- residualizeWithCovariates(reusable_genotype$covariate_adjustment,
                                                    gene_vector$phenotype_value_std)
    }
    names(expression_vector) <- sample_ids

    phenotype_chr <- normalize_chr(gene_meta$chromosome)
    phenotype_start <- gene_meta$phenotype_pos - cis_distance
    phenotype_end <- gene_meta$phenotype_pos + cis_distance
    variant_idx <- reusable_genotype$variant_info$chromosome == phenotype_chr &
      reusable_genotype$variant_info$position >= phenotype_start &
      reusable_genotype$variant_info$position <= phenotype_end
    variant_idx[is.na(variant_idx)] <- FALSE
    if (!any(variant_idx, na.rm = TRUE)) {
      message(
        "Reusable genotype matrix has no variants overlapping ",
        phenotype_chr,
        ":",
        phenotype_start,
        "-",
        phenotype_end,
        "; using per-phenotype genotype extraction"
      )
      message(
        "Reusable genotype chromosomes: ",
        paste(utils::head(unique(reusable_genotype$variant_info$chromosome), 10), collapse = ", ")
      )
    } else {
      gt_hat <- reusable_genotype$gt_hat[, variant_idx, drop = FALSE]
      variant_names <- reusable_genotype$variant_names[variant_idx]
      reportVariantDirectionCounts(
        variant_names,
        gene_meta,
        context = "reusable genotype matrix subset",
        output_file = variant_report_file,
        requested_window_start = phenotype_start,
        requested_window_end = phenotype_end
      )
      explained_variants <- reusable_genotype$covariate_explained_variants %>%
        dplyr::filter(variant_id %in% variant_names)
      fitted <- timed_fit(gt_hat, expression_vector, legacy = reusable_genotype$dense,
                              covariate_explained_variants = explained_variants)
      fitted$variant_id = variant_names
      fitted$cohort_maf <- reusable_genotype$cohort_maf[variant_idx]
      return(fitted)
    }
  }

  # Align covariates to phenotype sample order and add an intercept column.
  covariates_matrix = cbind(covariates[gene_vector$genotype_id,], 1)
  phenotype_start <- gene_meta$phenotype_pos - cis_distance
  phenotype_end <- gene_meta$phenotype_pos + cis_distance
  
  if (is.character(genotype_file)) {
  #Import genotype matrix
  message('Loading Genotype Matrix')
  genotype_matrix = eQTLUtils::extractGenotypeMatrixFromDosage(
    chr = gene_meta$chromosome, 
    start = phenotype_start,
    end = phenotype_end,
    dosage_file = genotype_file) %>% 
    filterMAF(ancestry_data,MAF_threshold = MAF,variant_list = variant_list) 
   } else {
    message('Using pre-loaded genotype matrix')
    genotype_matrix <- genotype_file
    }
  # Residualize expression with the same covariates and intercept.
  message('Residualizing gene expression')
  adjustment <- createCovariateAdjustment(covariates_matrix)
  expression_values <- gene_vector$phenotype_value_std
  expression_vector <- residualizeWithCovariates(adjustment, expression_values)
  names(expression_vector) = gene_vector$genotype_id

  if (!is.null(additional_genotypes)) {
   message('Merging with provided additional genotypes')
   genotype_matrix <- MergeAdditionalGenotypes(genotype_matrix,additional_genotypes)
   message('Subsetting to samples in the additional genotype file')
   merged_samples <- colnames(genotype_matrix)
   intersecting_samples <- intersect(merged_samples,as.character(names(expression_vector)))

   expression_vector <-  expression_vector[intersecting_samples, , drop = FALSE]
   names(expression_vector) <- rownames(expression_vector)
   covariates_matrix = covariates_matrix[names(expression_vector),]
   adjustment <- createCovariateAdjustment(covariates_matrix)
   subset_gene_vector <- gene_vector %>% filter(genotype_id %in% intersecting_samples)
   expression_values <- subset_gene_vector$phenotype_value_std
   expression_vector <- residualizeWithCovariates(adjustment, expression_values)
   names(expression_vector) <- subset_gene_vector$genotype_id  
  }


  # Restrict genotype columns to the samples present in molecular data.
  message('Subsetting genotype matrix to individuals in the molecular data')
  gt_matrix = genotype_matrix[,names(expression_vector)]
  rm(genotype_matrix)
  gc()
 # #Exclude variants with no alternative alleles
  #message('Filtering variants with no alt alleles')
  #gt_matrix = gt_matrix[rowSums(round(gt_matrix,0), na.rm = TRUE) != 0,]
    
  #message('Replacing missing genotypes ')
  ##Replace missing values with row means
  #gt_matrix = t(gt_matrix) %>% zoo::na.aggregate() %>% t()

  ##Standardise genotypes
  #message('Standardizing genotypes')
  #gt_std = t(gt_matrix - apply(gt_matrix, 1, mean))
  #gt_hat = hat %*% gt_std
 # subset genotype matrix to individuals that are in expression_vector
  # Exclude variants with no alt alleles
  nonzero_idx <- which(rowSums(round(gt_matrix,0), na.rm = TRUE) != 0)
  if (length(nonzero_idx) == 0) stop("No variant with alt alleles")
  gt_matrix <- gt_matrix[nonzero_idx, , drop = FALSE]
  if (isTRUE(getOption("susier.memory_debug", FALSE))) {
    message("DEBUG: before imputation (local objects and sizes):")
    debug_list_locals(environment(), top_n = 100)
    debug_print_specific(environment(), var_names = c("gt_matrix", "covariates_matrix"))
    message("GC before imputation:")
    print(gc())
  }
  # Impute missing values using row means and keep variant IDs before matrix
  # operations drop or reorder row names.
  all_missing <- if (genotypes_imputed) rep(FALSE, nrow(gt_matrix)) else {
    matrixStats::rowCounts(gt_matrix, value = NA_real_) +
      matrixStats::rowCounts(gt_matrix, value = -1, na.rm = TRUE) == ncol(gt_matrix)
  }
  if (!genotypes_imputed) {
    gt_matrix <- if (is.null(genotype_imputer)) {
      impute_matrix_rowmean(gt_matrix, missing_value = -1)
    } else {
      genotype_imputer(gt_matrix, missing_value = -1)
    }
  }
  if (!nrow(gt_matrix)) {
    stop("No variant with alt alleles")
  }
  variant_names <- rownames(gt_matrix)
  cohort_maf <- cohortMafFromDosages(gt_matrix)
  # The native imputer removes entirely missing rows. The R imputer keeps
  # zero-filled rows, whose observed frequency must remain unknown.
  if (length(all_missing) == nrow(gt_matrix)) cohort_maf[all_missing] <- NA_real_
  reportVariantDirectionCounts(
    variant_names,
    gene_meta,
    context = "fine-mapping genotype matrix",
    output_file = variant_report_file,
    requested_window_start = phenotype_start,
    requested_window_end = phenotype_end
  )
  # Standardize genotypes (subtract row means)

  message('Standardizing genotypes')
  gt_std <- standardize_rows(gt_matrix)   # still variants x samples
  rm(gt_matrix)
  gc()
  # The fit uses samples by variants; apply the small factor in that orientation.
  if (isTRUE(getOption("susier.memory_debug", FALSE))) {
    message("DEBUG: before computing gt_hat (local objects and sizes):")
    debug_list_locals(environment(), top_n = 100)
    debug_print_specific(environment(), var_names = c("gt_std", "covariates_matrix", "adjustment"))
    # Estimate memory needed for gt_hat if gt_std exists
    if (exists("gt_std", envir = environment()) && !is.null(dim(get("gt_std", envir = environment())))) {
      p <- dim(get("gt_std", envir = environment()))[1]   # variants
      n <- dim(get("gt_std", envir = environment()))[2]   # samples
      est_bytes <- as.numeric(n) * as.numeric(p) * 8
      message(sprintf("Estimated gt_hat memory: %s bytes (%.2f MB)", format(est_bytes, big.mark = ","), est_bytes / 1024^2))
    }
    message("GC before gt_hat:")
    print(gc())
  }
  message('Computing gt_hat')
  gt_hat <- residualizeWithCovariates(adjustment, t(gt_std), variant_ids = variant_names)
  legacy <- !is.null(adjustment$state$hat)
  if (legacy) {
    # Genotypes can trigger the guard after the response was residualized.
    # Repeat the response with the same original dense matrix in that case.
    expression_names <- names(expression_vector)
    expression_vector <- residualizeWithCovariates(adjustment, expression_values)
    names(expression_vector) <- expression_names
  }
  explained_variants <- adjustment$state$covariate_explained_variants
  rm(gt_std)
  rm(adjustment)
  gc()


  # Fit the SuSiE model on residualized genotypes and expression.
  fitted <- timed_fit(gt_hat, expression_vector, legacy = legacy,
                          covariate_explained_variants = explained_variants)
  fitted$variant_id = variant_names
  fitted$cohort_maf <- cohort_maf
  return(fitted)
}

# Convert a susieR fitted object into output tables for credible sets, variants,
# and per-effect log-Bayes factors.
extractResults <- function(susie_object){
  credible_sets = susie_object$sets$cs
  phenotype_residual_sd <- if (is.null(susie_object$phenotype_residual_SD)) {
    NA_real_
  } else {
    as.double(susie_object$phenotype_residual_SD)
  }
  flagged_variant_ids <- if (is.null(susie_object$covariate_explained_variants)) {
    character()
  } else {
    susie_object$covariate_explained_variants$variant_id
  }
  cs_list = list()
  susie_object$sets$purity = dplyr::as_tibble(susie_object$sets$purity) %>%
    dplyr::mutate(
      cs_id = rownames(susie_object$sets$purity),
      cs_size = NA,
      cs_log10bf = NA,
      overlapped = NA,
      has_covariate_explained_variant = FALSE
    )
  added_variants = c()
  for (index in seq_along(credible_sets)){
    cs_variants = credible_sets[[index]]
    cs_id = susie_object$sets$cs_index[[index]]
    susie_object$sets$purity$has_covariate_explained_variant[index] <-
      any(susie_object$variant_id[cs_variants] %in% flagged_variant_ids)

    is_overlapped = any(cs_variants %in% added_variants)
    susie_object$sets$purity$overlapped[index] = is_overlapped
    susie_object$sets$purity$cs_size[index] = length(cs_variants)
    susie_object$sets$purity$cs_log10bf[index] = log10(exp(susie_object$lbf[cs_id]))
    if (!is_overlapped) {
      cs_list[[index]] = dplyr::tibble(cs_id = paste0("L", cs_id),
                                       .predictor_index = cs_variants)
      added_variants = append(added_variants, cs_variants)
    }
  }
  df = purrr::map_df(cs_list, identity)

  # Extract purity values for all non-overlapping credible sets.
  purity_res = susie_object$sets$purity

  #Sometimes all the PIP values are 0 and there are no purity values, then skip this step
  if(nrow(purity_res) > 0){
    purity_df = dplyr::as_tibble(purity_res) %>%
      dplyr::filter(!overlapped) %>%
      dplyr::mutate(
        cs_avg_r2 = mean.abs.corr^2,
        cs_min_r2 = min.abs.corr^2,
        low_purity = min.abs.corr < 0.5
      )  %>%
      dplyr::select(cs_id, cs_log10bf, cs_avg_r2, cs_min_r2, cs_size, low_purity,
                    has_covariate_explained_variant)
  } else{
    purity_df = dplyr::tibble()
  }
  
  # Extract posterior means, standard errors, and lbf variables.
  mean_vec = susieR::susie_get_posterior_mean(susie_object)
  sd_vec = susieR::susie_get_posterior_sd(susie_object)
  
  #Extract matrices
  alpha_mat = t(susie_object$alpha)
  colnames(alpha_mat) = paste0("alpha", seq_len(ncol(alpha_mat)))
  
  mu_mat = t(susie_object$mu)
  colnames(mu_mat) = paste0("mu_", seq_len(ncol(mu_mat)))
  
  mu2_mat = t(susie_object$mu2)
  colnames(mu2_mat) = paste0("mu2_", seq_len(ncol(mu2_mat)))
  
  lbf_variable_mat = t(susie_object$lbf_variable)
  colnames(lbf_variable_mat) = paste0("lbf_variable", seq_len(ncol(lbf_variable_mat)))
  # Preserve the native variant PIP, including the fit's own prior threshold.
  # A credible-set assignment does not restrict the effects used in this PIP.
  posterior_df = dplyr::tibble(variant_id = rownames(alpha_mat), 
                               .predictor_index = seq_len(nrow(alpha_mat)),
                               pip = susie_object$pip,
                               z = susie_object$z[,1],
                               posterior_mean = mean_vec, 
                               posterior_sd = sd_vec,
                               X_column_scale_factors = susie_object$X_column_scale_factors,
                               phenotype_residual_SD = phenotype_residual_sd) %>%
                 dplyr::bind_cols(purrr::map(list(alpha_mat, mu_mat, mu2_mat), dplyr::as_tibble))
  lbf_df = dplyr::tibble(variant_id = rownames(lbf_variable_mat)) %>%
    dplyr::bind_cols(dplyr::as_tibble(lbf_variable_mat))

  if(nrow(df) > 0 && nrow(purity_df) > 0 && ncol(lbf_df) > 10){ #ncol(lbf_df) <= 10 only if the number of variants in the region is < 10
    cs_df = purity_df
    # Variant labels can repeat across physical predictors. Native CS members
    # use predictor column indices, so keep these indices until after the join.
    variant_df = dplyr::left_join(posterior_df, df, by = ".predictor_index") %>%
      dplyr::select(-.predictor_index) %>%
      dplyr::left_join(cs_df, by = "cs_id") %>%
      dplyr::mutate(has_covariate_explained_variant =
                      dplyr::coalesce(has_covariate_explained_variant, FALSE))
  } else{
    cs_df = NULL
    variant_df = NULL
    lbf_df = NULL
  }

  return(list(cs_df = cs_df, variant_df = variant_df, lbf_df = lbf_df))
}

# Build connected components across credible sets and keep the highest-PIP
# representative variant per component.
make_connected_components_from_cs <- function(susie_all_df, z_threshold = 3, cs_size_threshold = 10) {
  # Filter the credible sets by Z-score and size
  susie_filt_all <- susie_all_df %>%
    dplyr::group_by(group_id) %>%
    dplyr::mutate(max_abs_z = max(abs(z))) %>%
    dplyr::filter(max_abs_z > z_threshold, cs_size < cs_size_threshold) %>%
    dplyr::ungroup()
  
  susie_highest_pip_per_cc <- data.frame()
  
  uniq_groups = susie_filt_all$group_id %>% base::unique()
  # make the ranges object in order to find overlaps
  for (uniq_group in uniq_groups) {
    susie_filt <- susie_filt_all %>% dplyr::filter(group_id == uniq_group)
    message("Processing CC of group_id: ", uniq_group)

    cs_ranges = GenomicRanges::GRanges(
      seqnames = susie_filt$chromosome,
      ranges = IRanges::IRanges(start = susie_filt$position, end = susie_filt$position),
      strand = "*",
      mcols = data.frame(molecular_trait_id = susie_filt$molecular_trait_id, variant_id = susie_filt$variant, group_id = susie_filt$group_id)
    )
    
    # find overlaps and remove the duplicated
    olaps <- GenomicRanges::findOverlaps(cs_ranges, cs_ranges, ignore.strand = TRUE) %>%
      GenomicRanges::as.data.frame() %>%
      dplyr::filter(queryHits <= subjectHits)
    
    # change variant sharing into credible set sharing 
    # not to have multiple connected components of variants but credible sets
    olaps <- olaps %>% dplyr::mutate(cs_mol_1 = cs_ranges$mcols.molecular_trait_id[queryHits], cs_mol_2 = cs_ranges$mcols.molecular_trait_id[subjectHits])
    edge_list <- olaps %>% dplyr::select(cs_mol_1, cs_mol_2) %>% BiocGenerics::unique() %>% base::as.matrix()
    
    # make the graph of connected components
    g <- igraph::graph_from_edgelist(edge_list, directed = FALSE)
    g_cc <- igraph::components(g)
    
    # turn connected components graph into data frame
    cc_df <- data.frame(cc_membership_no = g_cc$membership, 
                        molecular_trait_id = g_cc$membership %>% names()) 
    
    susie_highest_pip_per_cc_temp <- susie_filt %>% 
      dplyr::left_join(cc_df, by = "molecular_trait_id") %>% 
      dplyr::group_by(cc_membership_no) %>% 
      dplyr::arrange(-pip) %>% 
      dplyr::slice(1) %>% 
      dplyr::ungroup()
    
    susie_highest_pip_per_cc <- susie_highest_pip_per_cc %>% base::rbind(susie_highest_pip_per_cc_temp)
  }
  
  return(susie_highest_pip_per_cc)
}

# Assign n items into fixed-size batches.
splitIntoBatches <- function(n, batch_size){
  n_batches = ceiling(n/batch_size)
  batch_ids = rep(seq(1:n_batches), each = batch_size)[1:n]
  return(batch_ids)
}

# Select the items belonging to one chunk among n_chunks total chunks.
splitIntoChunks <- function(chunk_number, n_chunks, n_total){
  chunk_size = max(1,floor(n_total/(n_chunks)))
  batches = splitIntoBatches(n_total,chunk_size)
  batches[batches > n_chunks] = seq(from = 1, to = length(batches[batches > n_chunks]))
  selected_batch = batches == chunk_number
  return(selected_batch)
}
