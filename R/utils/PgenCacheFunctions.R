# Prepared genotype caches are local to one batch, sample order, and covariate
# design. Physical PGEN indices identify columns even when variant names repeat.
.pgenCacheSamples <- function(reader, sample_ids) {
  if (!is.environment(reader) || !inherits(reader, "PgenReader") || isTRUE(reader$closed)) {
    stop("Prepared cache requires an open PgenReader", call. = FALSE)
  }
  if (!is.character(sample_ids) || !length(sample_ids) || anyNA(sample_ids) ||
      any(!nzchar(sample_ids)) || anyDuplicated(sample_ids) ||
      !all(sample_ids %in% reader$sample_ids)) {
    stop("Prepared cache sample_ids must be unique IDs in the PgenReader", call. = FALSE)
  }
  invisible(TRUE)
}

planPgenGenotypeCache <- function(reader, worklist, sample_ids, mode = "auto",
                                  cache_memory_gb = 0, job_memory_gb = 16,
                                  eligible = TRUE, reason = NULL) {
  .pgenCacheSamples(reader, sample_ids)
  if (!is.character(mode) || length(mode) != 1L || is.na(mode) ||
      !mode %in% c("auto", "whole", "rolling", "none")) {
    stop("Genotype cache mode must be auto, whole, rolling, or none", call. = FALSE)
  }
  valid_memory <- function(value, positive = FALSE) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && if (positive) value > 0 else value >= 0
  }
  if (!valid_memory(cache_memory_gb) || !valid_memory(job_memory_gb, TRUE) ||
      cache_memory_gb > job_memory_gb) {
    stop("Invalid genotype cache or job memory allocation", call. = FALSE)
  }
  if (!is.logical(eligible) || length(eligible) != 1L || is.na(eligible)) {
    stop("Cache eligibility must be TRUE or FALSE", call. = FALSE)
  }
  required <- c("outcome_id", "chromosome", "window_start", "window_end")
  if (!is.data.frame(worklist) || !nrow(worklist) ||
      !all(required %in% names(worklist)) || anyNA(worklist$outcome_id) ||
      any(!nzchar(as.character(worklist$outcome_id))) || anyDuplicated(worklist$outcome_id)) {
    stop("Cache worklist requires unique outcomes and exact feature windows", call. = FALSE)
  }
  windows <- lapply(seq_len(nrow(worklist)), function(i) {
    .pgenWindowIndices(reader, worklist$chromosome[[i]],
      worklist$window_start[[i]], worklist$window_end[[i]], allow_empty = TRUE)
  })
  names(windows) <- as.character(worklist$outcome_id)
  union_indices <- sort(unique(as.integer(unlist(windows, use.names = FALSE))))
  requested <- sum(lengths(windows))
  max_window <- max(lengths(windows))
  rolling_reads <- length(windows[[1L]])
  if (length(windows) > 1L) {
    for (i in seq.int(2L, length(windows))) {
      rolling_reads <- rolling_reads + length(setdiff(windows[[i]], windows[[i - 1L]]))
    }
  }
  whole_saved <- if (requested) 1 - length(union_indices) / requested else 0
  rolling_saved <- if (requested) 1 - rolling_reads / requested else 0
  block_size <- min(256L, max_window)
  # Four small input/output buffers cover block reading, centering, transpose,
  # and covariate adjustment. The fit's working matrices are separate.
  buffer_bytes <- 4 * 8 * length(sample_ids) * block_size
  whole_payload <- 8 * length(sample_ids) * length(union_indices)
  rolling_payload <- 8 * length(sample_ids) * max_window
  whole_bytes <- whole_payload + buffer_bytes
  rolling_bytes <- rolling_payload + buffer_bytes
  budget <- if (cache_memory_gb > 0) cache_memory_gb * 1e9 else .2 * job_memory_gb * 1e9
  fits <- function(bytes) bytes <= budget + max(1, budget) * 8 * .Machine$double.eps
  selected_mode <- mode
  selection_reason <- "explicit cache mode"
  if (mode == "none") {
    selection_reason <- "caching disabled"
  } else if (!eligible) {
    selection_reason <- if (is.null(reason)) "input preparation is not cache compatible" else as.character(reason)[[1L]]
    if (mode != "auto") stop("Requested genotype cache is unavailable: ", selection_reason, call. = FALSE)
    selected_mode <- "none"
  } else if (!length(union_indices)) {
    selection_reason <- "no variants in planned windows"
    if (mode != "auto") stop(selection_reason, call. = FALSE)
    selected_mode <- "none"
  } else if (mode == "auto") {
    if (fits(whole_bytes) && whole_saved >= .5) {
      selected_mode <- "whole"
      selection_reason <- "whole batch fits the budget and avoids at least 50% of reads"
    } else if (fits(rolling_bytes) && rolling_saved >= .5) {
      selected_mode <- "rolling"
      selection_reason <- "rolling window fits the budget and avoids at least 50% of reads"
    } else {
      selected_mode <- "none"
      selection_reason <- "cache memory or repeated-read savings do not meet the automatic rules"
    }
  } else if (!fits(if (mode == "whole") whole_bytes else rolling_bytes)) {
    stop("Requested genotype cache exceeds the cache memory budget", call. = FALSE)
  }
  list(mode = selected_mode, requested_mode = mode, reason = selection_reason,
    pgen_path = reader$pgen_path,
    budget_bytes = budget, sample_ids = sample_ids, sample_count = length(sample_ids),
    union_indices = union_indices, window_indices = windows,
    requested_reads = requested, readtotal = requested,
    union_count = length(union_indices), max_window = max_window,
    rolling_reads = rolling_reads, whole_saved_fraction = whole_saved,
    rolling_saved_fraction = rolling_saved, block_size = block_size,
    whole_payload_bytes = whole_payload, rolling_payload_bytes = rolling_payload,
    whole_bytes = whole_bytes, rolling_bytes = rolling_bytes,
    preparation_buffer_bytes = buffer_bytes,
    estimated_bytes = switch(selected_mode, whole = whole_bytes,
                             rolling = rolling_bytes, none = 0))
}

createPgenGenotypeCache <- function(reader, plan, sample_ids, covariates_matrix) {
  .pgenCacheSamples(reader, sample_ids)
  if (!is.list(plan) || !is.character(plan$mode) || length(plan$mode) != 1L ||
      !plan$mode %in% c("whole", "rolling") ||
      !identical(sample_ids, plan$sample_ids)) {
    stop("Prepared cache requires a whole or rolling plan with matching samples", call. = FALSE)
  }
  if (!identical(plan$pgen_path, reader$pgen_path)) {
    stop("Prepared cache PGEN input does not match the cache plan", call. = FALSE)
  }
  if (!is.matrix(covariates_matrix) || !is.numeric(covariates_matrix) ||
      nrow(covariates_matrix) != length(sample_ids) || !ncol(covariates_matrix) ||
      !identical(rownames(covariates_matrix), sample_ids) ||
      any(!is.finite(covariates_matrix))) {
    stop("Prepared cache covariate design must match sample order and contain finite values", call. = FALSE)
  }
  cache <- new.env(parent = emptyenv())
  class(cache) <- "PgenGenotypeCache"
  cache$reader <- reader
  cache$pgen_path <- reader$pgen_path
  cache$plan <- plan
  cache$mode <- plan$mode
  cache$sample_ids <- sample_ids
  cache$covariates_matrix <- covariates_matrix
  cache$projection <- covariates_matrix %*% solve(crossprod(covariates_matrix))
  cache$blocks <- list()
  cache$excluded_indices <- integer()
  cache$initialized <- FALSE
  cache$cleared <- FALSE
  cache$variants_read <- 0
  cache$variants_prepared <- 0
  cache$cache_hits <- 0
  cache$exclusion_hits <- 0
  cache$variants_evicted <- 0
  cache$guard_recomputations <- 0
  cache$peak_resident_bytes <- 0
  cache
}

.pgenCacheEmptyDiagnostics <- function() {
  tibble::tibble(variant_id = character(), input_scale = double(),
    residual_sd = double(), residual_sd_ratio = double())
}

.pgenCacheResidentIndices <- function(cache) {
  as.integer(unlist(lapply(cache$blocks, `[[`, "indices"), use.names = FALSE))
}

.pgenCachePrepare <- function(cache, indices) {
  if (!length(indices)) return(invisible(NULL))
  starts <- seq.int(1L, length(indices), by = cache$plan$block_size)
  for (start in starts) {
    requested <- indices[seq.int(start, min(length(indices), start + cache$plan$block_size - 1L))]
    decoded <- readPgenVariants(cache$reader, requested, sample_ids = cache$sample_ids,
      as_double = TRUE, meanimpute = TRUE)
    cache$variants_read <- cache$variants_read + length(requested)
    mat <- decoded$genotypes
    retained <- decoded$variant_indices
    if (length(retained)) {
      nonzero <- rowSums(round(mat, 0), na.rm = TRUE) != 0
      mat <- mat[nonzero, , drop = FALSE]
      retained <- retained[nonzero]
    }
    excluded <- setdiff(requested, retained)
    cache$excluded_indices <- union(cache$excluded_indices, excluded)
    if (!length(retained)) next
    values <- t(mat - rowMeans(mat))
    residuals <- values - cache$projection %*% crossprod(cache$covariates_matrix, values)
    input_scale <- pmax(matrixStats::colSds(values), abs(colMeans(values)))
    residual_sd <- matrixStats::colSds(residuals)
    flagged <- which(input_scale > 0 & residual_sd <= 1e-8 * input_scale)
    diagnostics <- tibble::tibble(
      variant_id = cache$reader$variants$variant_id[retained[flagged]],
      input_scale = unname(input_scale[flagged]),
      residual_sd = unname(residual_sd[flagged]),
      residual_sd_ratio = unname(residual_sd[flagged] / input_scale[flagged]))
    cache$blocks[[length(cache$blocks) + 1L]] <- list(
      indices = as.integer(retained), gt_hat = residuals,
      flagged_indices = as.integer(retained[flagged]), diagnostics = diagnostics)
    cache$variants_prepared <- cache$variants_prepared + length(retained)
    resident_bytes <- 8 * length(cache$sample_ids) * length(.pgenCacheResidentIndices(cache))
    cache$peak_resident_bytes <- max(cache$peak_resident_bytes, resident_bytes)
  }
  invisible(NULL)
}

.pgenCacheEvict <- function(cache, required) {
  kept <- list()
  for (block in cache$blocks) {
    retain <- block$indices %in% required
    cache$variants_evicted <- cache$variants_evicted + sum(!retain)
    if (!any(retain)) next
    if (!all(retain)) {
      block$indices <- block$indices[retain]
      block$gt_hat <- block$gt_hat[, retain, drop = FALSE]
      keep_flags <- block$flagged_indices %in% block$indices
      block$flagged_indices <- block$flagged_indices[keep_flags]
      block$diagnostics <- block$diagnostics[keep_flags, , drop = FALSE]
    }
    kept[[length(kept) + 1L]] <- block
  }
  cache$blocks <- kept
  invisible(NULL)
}

getPgenGenotypeCacheWindow <- function(cache, outcome_id) {
  if (!is.environment(cache) || !inherits(cache, "PgenGenotypeCache") || isTRUE(cache$cleared)) {
    stop("Prepared genotype cache is missing or cleared", call. = FALSE)
  }
  .pgenCacheSamples(cache$reader, cache$sample_ids)
  if (!identical(cache$pgen_path, cache$reader$pgen_path)) {
    stop("Prepared cache PGEN identity has changed", call. = FALSE)
  }
  if (!is.character(outcome_id) || length(outcome_id) != 1L || is.na(outcome_id) ||
      !outcome_id %in% names(cache$plan$window_indices)) {
    stop("Outcome is absent from the prepared cache plan", call. = FALSE)
  }
  required <- cache$plan$window_indices[[outcome_id]]
  resident <- .pgenCacheResidentIndices(cache)
  cache$cache_hits <- cache$cache_hits + length(intersect(required, resident))
  cache$exclusion_hits <- cache$exclusion_hits + length(intersect(required, cache$excluded_indices))
  if (cache$mode == "whole" && !cache$initialized) {
    .pgenCachePrepare(cache, cache$plan$union_indices)
    cache$initialized <- TRUE
  } else if (cache$mode == "rolling") {
    evicted_before <- cache$variants_evicted
    .pgenCacheEvict(cache, required)
    # Collect the old blocks after the eviction frame has released them and
    # before incoming variants allocate new genotype and adjustment matrices.
    if (cache$variants_evicted > evicted_before) invisible(gc())
    missing <- setdiff(required, c(.pgenCacheResidentIndices(cache), cache$excluded_indices))
    .pgenCachePrepare(cache, missing)
    cache$initialized <- TRUE
  }
  indices <- required[required %in% .pgenCacheResidentIndices(cache)]
  names <- cache$reader$variants$variant_id[indices]
  result <- matrix(0, nrow = length(cache$sample_ids), ncol = length(indices),
    dimnames = list(cache$sample_ids, names))
  diagnostics <- list()
  guard <- FALSE
  for (block in cache$blocks) {
    in_window <- which(block$indices %in% indices)
    if (!length(in_window)) next
    target <- match(block$indices[in_window], indices)
    result[, target] <- block$gt_hat[, in_window, drop = FALSE]
    flags <- block$flagged_indices %in% indices
    if (any(flags)) {
      guard <- TRUE
      diagnostics[[length(diagnostics) + 1L]] <- block$diagnostics[flags, , drop = FALSE]
    }
  }
  explained <- if (length(diagnostics)) dplyr::bind_rows(diagnostics) else .pgenCacheEmptyDiagnostics()
  # Diagnostic order follows physical PGEN order, matching the feature matrix.
  if (nrow(explained)) {
    explained <- explained[order(match(explained$variant_id, names)), , drop = FALSE]
  }
  list(gt_hat = result, variant_indices = as.integer(indices), variant_names = names,
    covariate_explained_variants = explained, guard_required = guard)
}

pgenGenotypeCacheSummary <- function(cache) {
  if (!is.environment(cache) || !inherits(cache, "PgenGenotypeCache")) {
    stop("Invalid prepared genotype cache", call. = FALSE)
  }
  resident <- length(.pgenCacheResidentIndices(cache))
  list(mode = cache$mode, variants_read = cache$variants_read,
    variants_prepared = cache$variants_prepared, cache_hits = cache$cache_hits,
    exclusion_hits = cache$exclusion_hits, variants_evicted = cache$variants_evicted,
    excluded_variants = length(cache$excluded_indices), resident_variants = resident,
    resident_bytes = 8 * length(cache$sample_ids) * resident,
    peak_resident_bytes = cache$peak_resident_bytes,
    guard_recomputations = cache$guard_recomputations)
}

clearPgenGenotypeCache <- function(cache) {
  if (!is.environment(cache) || !inherits(cache, "PgenGenotypeCache")) {
    stop("Invalid prepared genotype cache", call. = FALSE)
  }
  cache$blocks <- list()
  cache$cleared <- TRUE
  invisible(NULL)
}
