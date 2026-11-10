# Request memory from metadata. No genotype matrix is read or allocated here.
# All sizes use decimal GB, as in the WDL task memory input.
validatePgenBatchMemoryOptions <- function(cache_mode = "auto", cache_memory_gb = 0,
                                           memory_gb = 0, cache_eligible = TRUE) {
  if (!is.character(cache_mode) || length(cache_mode) != 1L || is.na(cache_mode) ||
      !cache_mode %in% c("auto", "whole", "rolling", "none")) {
    stop("Genotype cache mode must be auto, whole, rolling, or none", call. = FALSE)
  }
  valid_memory <- function(value) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value >= 0
  }
  if (!valid_memory(cache_memory_gb) || !valid_memory(memory_gb) ||
      memory_gb != floor(memory_gb) || memory_gb > .Machine$integer.max ||
      cache_memory_gb > .Machine$integer.max ||
      (memory_gb > 0 && cache_memory_gb > memory_gb)) {
    stop("Invalid cache or task memory: task memory must be a nonnegative integer; ",
      "cache memory must fit within a positive task override", call. = FALSE)
  }
  if (!is.logical(cache_eligible) || length(cache_eligible) != 1L || is.na(cache_eligible)) {
    stop("Cache eligibility must be TRUE or FALSE", call. = FALSE)
  }
  if (!cache_eligible && cache_mode %in% c("whole", "rolling")) {
    stop("Requested genotype cache is unavailable for these inputs", call. = FALSE)
  }
  invisible(TRUE)
}

.pgenMemoryCount <- function(value, label, minimum = 0) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < minimum || value != floor(value) ||
      value > .Machine$integer.max) {
    stop(label, " must be a finite whole count within the R/WDL integer range", call. = FALSE)
  }
  as.double(value)
}

.pgenMemoryIntegerGB <- function(value) {
  rounded <- ceiling(value)
  if (!is.finite(rounded) || rounded < 1 || rounded > .Machine$integer.max) {
    stop("Estimated task memory exceeds the R/WDL integer range", call. = FALSE)
  }
  as.integer(rounded)
}

# The fit allowance includes five copies of the largest double window and
# eight L-by-P arrays for each retained outcome (L = 10). Two sample-by-sample
# matrices cover the compatibility guard. Reserve_fraction solves the cache
# budget and safety allowance together, without iteration.
.pgenMemoryAllocation <- function(sample_count, max_window_variants, cache_bytes = 0,
                                    reserve_fraction = 0, minimum_memory_gb = 0,
                                    max_gene_requested_variants = 0,
                                    total_requested_variants = 0) {
  n <- .pgenMemoryCount(sample_count, "sample_count", 1)
  p <- .pgenMemoryCount(max_window_variants, "max_window_variants")
  gene_variants <- .pgenMemoryCount(max_gene_requested_variants, "max_gene_requested_variants")
  total_variants <- .pgenMemoryCount(total_requested_variants, "total_requested_variants")
  if (!is.numeric(cache_bytes) || length(cache_bytes) != 1L ||
      !is.finite(cache_bytes) || cache_bytes < 0 ||
      !is.numeric(reserve_fraction) || length(reserve_fraction) != 1L ||
      !is.finite(reserve_fraction) || reserve_fraction < 0 || reserve_fraction >= .8 ||
      !is.numeric(minimum_memory_gb) || length(minimum_memory_gb) != 1L ||
      !is.finite(minimum_memory_gb) || minimum_memory_gb < 0) {
    stop("Invalid cache memory allowance", call. = FALSE)
  }
  matrix_bytes <- 8 * n * p
  retained_fit_bytes <- 64 * 10 * gene_variants
  fit_bytes <- max(8e9, 5 * matrix_bytes + retained_fit_bytes)
  dense_guard_bytes <- 2 * 8 * n^2
  # Compact tables remain live until batch aggregation. This per-row allowance
  # is a planning proxy, not a bound on R object size.
  summary_bytes <- 512 * total_variants
  overhead_bytes <- max(2e9, summary_bytes)
  safety_fraction <- .25
  base_bytes <- fit_bytes + dense_guard_bytes + overhead_bytes
  request <- max(1.25 * (base_bytes + cache_bytes) / 1e9,
    1.25 * base_bytes / (1 - 1.25 * reserve_fraction) / 1e9, minimum_memory_gb)
  memory <- .pgenMemoryIntegerGB(request)
  list(matrix_bytes = matrix_bytes, fit_bytes = fit_bytes,
    retained_fit_bytes = retained_fit_bytes, summary_bytes = summary_bytes,
    dense_guard_bytes = dense_guard_bytes, overhead_bytes = overhead_bytes,
    safety_fraction = safety_fraction,
    cache_bytes = max(cache_bytes, reserve_fraction * memory * 1e9),
    estimated_memory_gb = memory)
}

planPgenBatchMemory <- function(reader, worklist, cache_mode = "auto",
                                 cache_memory_gb = 0, memory_gb = 0,
                                 cache_eligible = TRUE) {
  validatePgenBatchMemoryOptions(cache_mode, cache_memory_gb, memory_gb, cache_eligible)
  .pgenCacheSamples(reader, reader$sample_ids)
  sample_count <- .pgenMemoryCount(reader$raw_sample_count, "sample_count", 1)
  raw_ids <- reader$raw_sample_ids
  if (!is.character(raw_ids) || length(raw_ids) != sample_count || anyNA(raw_ids) ||
      any(!nzchar(raw_ids)) || anyDuplicated(raw_ids) ||
      !all(reader$sample_ids %in% raw_ids)) {
    stop("PSAM sample count or ordered sample IDs are invalid", call. = FALSE)
  }
  if (!is.data.frame(worklist) || !nrow(worklist) || !"gene_id" %in% names(worklist) ||
      anyNA(worklist$gene_id) || any(!nzchar(as.character(worklist$gene_id)))) {
    stop("Memory worklist requires a nonempty gene_id for each outcome", call. = FALSE)
  }
  # Match run_pgen_susie_batch.R: process all outcomes of a gene together, in
  # first-seen gene order. Rolling reuse depends on this execution order.
  execution <- worklist[order(match(worklist$gene_id, unique(worklist$gene_id))), , drop = FALSE]
  # This view contains metadata only. Use the full PSAM cohort as an upper
  # bound even when the opened reader was configured for fewer samples.
  metadata <- new.env(parent = emptyenv())
  class(metadata) <- "PgenReader"
  metadata$closed <- FALSE
  metadata$pgen_path <- reader$pgen_path
  metadata$sample_ids <- raw_ids
  metadata$variants <- reader$variants
  metadata$position_lookup <- reader$position_lookup
  replan <- function(mode, job_gb, cap_gb = cache_memory_gb) {
    planPgenGenotypeCache(metadata, execution, raw_ids, mode = mode,
      cache_memory_gb = cap_gb, job_memory_gb = job_gb, eligible = cache_eligible)
  }
  # Count each physical index. Duplicate positions and IDs remain distinct.
  counts <- replan("none", max(1, cache_memory_gb))
  max_window <- .pgenMemoryCount(counts$max_window, "max_window_variants")
  union_count <- .pgenMemoryCount(counts$union_count, "union_variants")
  per_gene_windows <- split(lengths(counts$window_indices), as.character(execution$gene_id))
  max_gene_variants <- .pgenMemoryCount(max(vapply(per_gene_windows, sum, numeric(1))),
    "max_gene_requested_variants")
  max_gene_outcomes <- .pgenMemoryCount(max(lengths(per_gene_windows)), "max_gene_outcomes", 1)
  total_variants <- .pgenMemoryCount(sum(lengths(counts$window_indices)),
    "total_requested_variants")
  if (union_count < max_window || union_count > nrow(reader$variants)) {
    stop("Variant counts are inconsistent with PVAR metadata", call. = FALSE)
  }
  reuse <- max(counts$whole_saved_fraction, counts$rolling_saved_fraction) >= .5
  reserve_cache <- cache_eligible && union_count > 0 &&
    (cache_mode %in% c("whole", "rolling") || (cache_mode == "auto" && reuse))
  candidate_bytes <- if (!reserve_cache) 0 else if (cache_mode == "rolling") {
    counts$rolling_bytes
  } else counts$whole_bytes
  reserve_fraction <- if (reserve_cache && cache_memory_gb == 0 &&
    cache_mode %in% c("auto", "rolling")) .2 else 0
  minimum_gb <- cache_memory_gb
  if (reserve_cache && cache_memory_gb == 0 && cache_mode %in% c("whole", "rolling")) {
    minimum_gb <- max(minimum_gb, 5 * candidate_bytes / 1e9)
  }
  allowance <- if (reserve_cache && cache_memory_gb > 0) cache_memory_gb * 1e9 else candidate_bytes
  # Start with the union candidate, then replan with the resulting budget. At
  # most three passes can move whole -> rolling -> none; no pass upgrades it.
  estimate <- .pgenMemoryAllocation(sample_count, max_window, allowance,
    reserve_fraction, minimum_gb, max_gene_variants, total_variants)
  for (pass in seq_len(3L)) {
    predicted <- replan(cache_mode, estimate$estimated_memory_gb)
    allowance <- if (reserve_cache && cache_memory_gb > 0) {
      cache_memory_gb * 1e9
    } else predicted$estimated_bytes
    updated <- .pgenMemoryAllocation(sample_count, max_window, allowance,
      reserve_fraction, minimum_gb, max_gene_variants, total_variants)
    if (identical(updated$estimated_memory_gb, estimate$estimated_memory_gb)) {
      estimate <- updated
      break
    }
    estimate <- updated
  }
  actual_memory <- if (memory_gb > 0) as.integer(memory_gb) else estimate$estimated_memory_gb
  predicted <- replan(cache_mode, actual_memory)
  actual_cache <- if (reserve_cache && cache_memory_gb > 0) {
    cache_memory_gb * 1e9
  } else max(predicted$estimated_bytes, reserve_fraction * actual_memory * 1e9)
  list(memory_gb = actual_memory, cache_mode = predicted$mode,
    model_version = "pgen-batch-memory-v1", sample_count = as.integer(sample_count),
    max_window_variants = as.integer(max_window), union_variants = as.integer(union_count),
    max_gene_requested_variants = as.integer(max_gene_variants),
    max_gene_outcomes = as.integer(max_gene_outcomes),
    total_requested_variants = as.integer(total_variants),
    matrix_bytes = estimate$matrix_bytes, fit_bytes = estimate$fit_bytes,
    retained_fit_bytes = estimate$retained_fit_bytes, summary_bytes = estimate$summary_bytes,
    dense_guard_bytes = estimate$dense_guard_bytes,
    cache_bytes = actual_cache, cache_payload_bytes = predicted$estimated_bytes,
    cache_budget_bytes = predicted$budget_bytes, overhead_bytes = estimate$overhead_bytes,
    safety_fraction = estimate$safety_fraction,
    estimated_memory_gb = estimate$estimated_memory_gb)
}
