# Metadata and count tests must not decode genotypes or allocate fit matrices.
run_pgen_memory_tests <- function(repo_dir) {
  source(file.path(repo_dir, "tests", "test_helpers.R"), local = TRUE)
  functions_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH", file.path(repo_dir, "R", "utils"))
  for (filename in c("PgenFunctions.R", "PgenCacheFunctions.R", "PgenBatchFunctions.R",
                     "PgenMemoryFunctions.R")) {
    path <- file.path(functions_path, filename)
    if (file.exists(path)) source(path, local = TRUE)
  }
  run_named_test("memory options are validated before genotype preparation starts", {
    stopifnot(exists("validatePgenBatchMemoryOptions", mode = "function"))
    validatePgenBatchMemoryOptions()
    validatePgenBatchMemoryOptions(cache_memory_gb = .25, memory_gb = 1)
    expect_error_message(validatePgenBatchMemoryOptions(memory_gb = .5), "memory")
    expect_error_message(validatePgenBatchMemoryOptions(cache_memory_gb = 2,
      memory_gb = 1), "memory")
    expect_error_message(validatePgenBatchMemoryOptions(cache_mode = "whole",
      cache_eligible = FALSE), "unavailable")
  })
  run_named_test("fit memory keeps the eight GB floor and adds guard and safety", {
    stopifnot(exists(".pgenMemoryAllocation", mode = "function"))
    small <- .pgenMemoryAllocation(8, 5)
    expect_identical_value(small$estimated_memory_gb, 13L)
    stopifnot(small$matrix_bytes == 320, small$fit_bytes == 8e9,
      small$dense_guard_bytes == 1024, small$overhead_bytes == 2e9,
      small$safety_fraction == .25)
    ordinary <- .pgenMemoryAllocation(10000, 20000)
    expect_identical_value(ordinary$estimated_memory_gb, 15L)
    large <- .pgenMemoryAllocation(20000, 30000)
    expect_identical_value(large$estimated_memory_gb, 41L)
    stopifnot(large$matrix_bytes == 4.8e9, large$fit_bytes == 24e9,
      large$dense_guard_bytes == 6.4e9)
    reserved <- .pgenMemoryAllocation(20000, 30000, reserve_fraction = .2)
    expect_identical_value(reserved$estimated_memory_gb, 54L)
    stopifnot(reserved$cache_bytes >= .2 * 54e9)
  })
  run_named_test("count formula rejects invalid and unrepresentable requests", {
    for (bad in c(-1, .5, NA_real_, Inf, NaN, .Machine$integer.max + 1)) {
      expect_error_message(.pgenMemoryAllocation(bad, 5), "sample_count")
      expect_error_message(.pgenMemoryAllocation(8, bad), "max_window_variants")
      expect_error_message(.pgenMemoryAllocation(8, 5,
        max_gene_requested_variants = bad), "max_gene_requested_variants")
      expect_error_message(.pgenMemoryAllocation(8, 5,
        total_requested_variants = bad), "total_requested_variants")
    }
    expect_error_message(.pgenMemoryAllocation(0, 5), "sample_count")
    expect_error_message(.pgenMemoryAllocation(.Machine$integer.max,
      .Machine$integer.max), "integer range")
  })
  run_named_test("retained outcomes add fit arrays without adding full input matrices", {
    one <- .pgenMemoryAllocation(10000, 20000,
      max_gene_requested_variants = 20000, total_requested_variants = 20000)
    repeated <- .pgenMemoryAllocation(10000, 20000,
      max_gene_requested_variants = 2000000, total_requested_variants = 2000000)
    expect_identical_value(one$estimated_memory_gb, 15L)
    expect_identical_value(repeated$estimated_memory_gb, 17L)
    stopifnot(one$retained_fit_bytes == 12800000,
      repeated$retained_fit_bytes == 1280000000, repeated$fit_bytes == 9.28e9,
      repeated$matrix_bytes == one$matrix_bytes)
  })
  run_named_test("batch summary allowance increases overhead for many requested variants", {
    result <- .pgenMemoryAllocation(8, 5, total_requested_variants = 10000000)
    expect_identical_value(result$estimated_memory_gb, 17L)
    stopifnot(result$summary_bytes == 5.12e9, result$overhead_bytes == 5.12e9)
  })
  directory <- tempfile("pgen-memory-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  samples <- paste0("s", seq_len(8L))
  files <- file.path(directory, c("calls.pgen", "calls.pvar", "calls.psam"))
  writer <- pgenlibr::NewPgenWriter(files[[1L]], 8L, 12L, nonref_flags = FALSE)
  tryCatch({
    for (i in seq_len(12L)) pgenlibr::AppendBiallelic(writer, c(0, 1, 2, 0, 1, 2, 1, 0))
  }, finally = pgenlibr::ClosePgenWriter(writer))
  positions <- c(10L, 100L, 140L, 110L, 120L, 130L, 150L, 160L, 170L, 180L, 190L, 190L)
  writeLines(c("#CHROM\tPOS\tID\tREF\tALT", paste(c(2L, rep(1L, 11L)),
    positions, "duplicate-id", "A", "G", sep = "\t")), files[[2L]])
  writeLines(c("#IID", samples), files[[3L]])
  reader <- do.call(openPgenReader, c(as.list(files), list(sample_ids = samples[1:6])))
  on.exit(closePgenReader(reader), add = TRUE)
  worklist <- data.frame(gene_id = c("a", "b", "a", "b"),
    outcome_id = LETTERS[1:4], chromosome = "1", window_start = c(100L, 130L, 110L, 120L),
    window_end = c(140L, 170L, 150L, 160L))
  # Any genotype decode is a failure. Metadata opening above is allowed.
  for (name in c("ReadList", "ReadIntList")) {
    suppressMessages(trace(name, tracer = quote(stop("Memory planning decoded genotypes")),
      where = asNamespace("pgenlibr"), print = FALSE))
  }
  on.exit({
    for (name in c("ReadList", "ReadIntList")) {
      suppressMessages(untrace(name, where = asNamespace("pgenlibr")))
    }
  }, add = TRUE)
  run_named_test("metadata planning uses full PSAM counts and physical variants", {
    result <- planPgenBatchMemory(reader, worklist)
    expect_identical_value(result$memory_gb, 17L)
    expect_identical_value(result$estimated_memory_gb, 17L)
    expect_identical_value(result$cache_mode, "whole")
    stopifnot(result$sample_count == 8, result$max_window_variants == 5,
      result$union_variants == 8, result$matrix_bytes == 320,
      result$cache_payload_bytes == 1792, result$cache_bytes >= .2 * 17e9,
      result$max_gene_requested_variants == 10, result$max_gene_outcomes == 2,
      result$total_requested_variants == 20, result$retained_fit_bytes == 6400,
      result$summary_bytes == 10240,
      is.character(result$model_version), nzchar(result$model_version))
    duplicates <- worklist[1:2, ]
    duplicates$window_start <- duplicates$window_end <- 190L
    duplicate_result <- planPgenBatchMemory(reader, duplicates, cache_mode = "none")
    stopifnot(duplicate_result$max_window_variants == 2,
      duplicate_result$union_variants == 2)
  })
  run_named_test("auto avoids cache inflation when windows have little reuse", {
    disjoint <- worklist[1:2, ]
    disjoint$window_start <- disjoint$window_end <- c(100L, 190L)
    result <- planPgenBatchMemory(reader, disjoint)
    expect_identical_value(result$cache_mode, "none")
    expect_identical_value(result$memory_gb, 13L)
    stopifnot(result$cache_bytes == 0, result$cache_payload_bytes == 0)
    unavailable <- planPgenBatchMemory(reader, worklist, cache_eligible = FALSE)
    expect_identical_value(unavailable$cache_mode, "none")
    expect_identical_value(unavailable$memory_gb, 13L)
    expect_error_message(planPgenBatchMemory(reader, worklist,
      cache_mode = "whole", cache_eligible = FALSE), "unavailable")
  })
  run_named_test("explicit budgets select rolling or none and protect forced modes", {
    rolling <- planPgenBatchMemory(reader, worklist, cache_memory_gb = 1696 / 1e9)
    expect_identical_value(rolling$cache_mode, "rolling")
    stopifnot(rolling$cache_payload_bytes == 1600, rolling$cache_bytes == 1696)
    none <- planPgenBatchMemory(reader, worklist, cache_memory_gb = 800 / 1e9)
    expect_identical_value(none$cache_mode, "none")
    expect_error_message(planPgenBatchMemory(reader, worklist,
      cache_mode = "whole", cache_memory_gb = 1696 / 1e9), "budget")
    forced <- planPgenBatchMemory(reader, worklist, cache_mode = "rolling")
    expect_identical_value(forced$cache_mode, "rolling")
    stopifnot(forced$cache_bytes >= .2 * forced$memory_gb * 1e9,
      forced$cache_payload_bytes <= .2 * forced$memory_gb * 1e9)
    capped <- planPgenBatchMemory(reader, worklist, cache_memory_gb = 30)
    stopifnot(capped$memory_gb >= 30L, capped$cache_bytes == 30e9)
  })
  run_named_test("positive task memory overrides the request without hiding the estimate", {
    overridden <- planPgenBatchMemory(reader, worklist, memory_gb = 4)
    expect_identical_value(overridden$memory_gb, 4L)
    expect_identical_value(overridden$estimated_memory_gb, 17L)
    stopifnot(overridden$cache_bytes >= .2 * 4e9)
    expect_error_message(planPgenBatchMemory(reader, worklist,
      cache_memory_gb = 5, memory_gb = 4), "memory")
  })
  run_named_test("a large union downgrades to rolling without inflating task memory", {
    # Scale metadata only; the small real PGEN handle is never decoded.
    large <- list2env(as.list.environment(reader, all.names = TRUE), parent = emptyenv())
    class(large) <- "PgenReader"
    large$raw_sample_ids <- large$sample_ids <- paste0("n", seq_len(20000L))
    large$raw_sample_count <- 20000L
    large$variants <- data.frame(pgen_index = seq_len(300000L), chromosome = "1",
      position = seq_len(300000L), variant_id = "duplicate")
    large$position_lookup <- list("1" = list(positions = seq_len(300000L),
      indices = seq_len(300000L)))
    group <- rep(seq_len(300L), each = 3L)
    many <- data.frame(gene_id = paste0("g", group),
      outcome_id = paste0("o", seq_along(group)), chromosome = "1",
      window_start = (group - 1L) * 1000L + 1L, window_end = group * 1000L)
    result <- planPgenBatchMemory(large, many)
    expect_identical_value(result$cache_mode, "rolling")
    expect_identical_value(result$memory_gb, 28L)
    stopifnot(result$union_variants == 300000L, result$max_window_variants == 1000L,
      result$cache_payload_bytes == 323840000, abs(result$cache_bytes - 5.6e9) < 1)
    forced <- planPgenBatchMemory(large, many, cache_mode = "whole")
    expect_identical_value(forced$cache_mode, "whole")
    stopifnot(forced$cache_payload_bytes <= .2 * forced$memory_gb * 1e9)
    expect_error_message(planPgenBatchMemory(large, many,
      cache_mode = "whole", memory_gb = 28), "budget")
  })
  run_named_test("bad inputs fail before a memory request is returned", {
    for (bad in c(-1, .5, NA_real_, Inf, NaN, .Machine$integer.max + 1)) {
      expect_error_message(planPgenBatchMemory(reader, worklist, memory_gb = bad), "memory")
    }
    expect_error_message(planPgenBatchMemory(reader, worklist, cache_mode = "bad"), "mode")
    expect_error_message(planPgenBatchMemory(reader, worklist,
      cache_memory_gb = -1), "memory")
    expect_error_message(planPgenBatchMemory(reader, worklist,
      cache_eligible = NA), "eligibility")
    invalid <- worklist
    invalid$gene_id[[1L]] <- NA_character_
    expect_error_message(planPgenBatchMemory(reader, invalid), "gene_id")
    invalid <- worklist
    invalid$window_start[[1L]] <- 0
    expect_error_message(planPgenBatchMemory(reader, invalid), "bounds")
    original_count <- reader$raw_sample_count
    reader$raw_sample_count <- 9L
    tryCatch(expect_error_message(planPgenBatchMemory(reader, worklist), "PSAM"),
      finally = reader$raw_sample_count <- original_count)
  })
  message("All metadata PGEN memory tests passed.")
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  benchmark_library <- Sys.getenv("BENCH_LIBRARY", "")
  if (nzchar(benchmark_library)) .libPaths(c(benchmark_library, .libPaths()))
  args <- commandArgs(trailingOnly = TRUE)
  repo_dir <- normalizePath(if (length(args)) args[[1L]] else getwd(), mustWork = TRUE)
  run_pgen_memory_tests(repo_dir)
}
