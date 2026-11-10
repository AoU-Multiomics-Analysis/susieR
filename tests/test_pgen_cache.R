# Real PGEN files test cache selection, physical variant mapping, and reuse.
run_pgen_cache_tests <- function(repo_dir) {
  source(file.path(repo_dir, "tests", "test_helpers.R"), local = TRUE)
  functions_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH", file.path(repo_dir, "R", "utils"))
  for (filename in c("SusieFunctions.R", "PgenFunctions.R", "PgenCacheFunctions.R")) {
    path <- file.path(functions_path, filename)
    if (file.exists(path)) source(path, local = TRUE)
  }
  directory <- tempfile("prepared-pgen-cache-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  samples <- paste0("s", seq_len(8L))
  selected <- samples[c(7L, 3L, 1L, 5L, 4L, 8L)]
  set.seed(721L)
  calls <- matrix(as.double(rbinom(12L * 8L, 2L, .35)), nrow = 12L)
  calls[2L, ] <- c(0, 0, 1, 1, 2, 2, 1, 0)
  calls[4L, 7L] <- NA_real_
  calls[5L, match(selected, samples)] <- NA_real_
  calls[6L, ] <- 0
  colnames(calls) <- samples
  positions <- c(10L, 100L, 140L, 110L, 120L, 130L, 150L, 160L, 170L, 180L, 190L, 200L)
  files <- file.path(directory, c("calls.pgen", "calls.pvar", "calls.psam"))
  writer <- pgenlibr::NewPgenWriter(files[[1L]], 8L, 12L, nonref_flags = FALSE)
  tryCatch({
    for (i in seq_len(nrow(calls))) pgenlibr::AppendBiallelic(writer, calls[i, ])
  }, finally = pgenlibr::ClosePgenWriter(writer))
  writeLines(c("#CHROM\tPOS\tID\tREF\tALT", paste(c(2L, rep(1L, 11L)),
    positions, paste0("v", seq_len(12L)), "A", "G", sep = "\t")), files[[2L]])
  writeLines(c("#IID", samples), files[[3L]])
  reader <- do.call(openPgenReader, as.list(files))
  on.exit(closePgenReader(reader), add = TRUE)
  worklist <- data.frame(outcome_id = LETTERS[1:4], chromosome = "1",
    window_start = c(100L, 110L, 120L, 130L),
    window_end = c(140L, 150L, 160L, 170L))
  design <- cbind(pc = c(-2, -1, 0, 1, 2, 3), intercept = 1)
  rownames(design) <- selected
  observe_reads <- function(code) {
    observer <- new.env(parent = emptyenv())
    observer$indices <- list()
    tracer <- substitute({
      .cache_read_observer <- captured
      .cache_read_observer$indices[[length(.cache_read_observer$indices) + 1L]] <- variant_subset
    }, list(captured = observer))
    suppressMessages(trace("ReadList", tracer = tracer,
      where = asNamespace("pgenlibr"), print = FALSE))
    on.exit(suppressMessages(untrace("ReadList", where = asNamespace("pgenlibr"))), add = TRUE)
    value <- force(code)
    list(value = value, indices = observer$indices)
  }
  run_named_test("automatic cache selection uses overlap and exact memory limits", {
    whole <- planPgenGenotypeCache(reader, worklist, selected)
    expect_identical_value(whole$mode, "whole")
    expect_identical_value(whole$window_indices$A, 2:6)
    expect_identical_value(whole$union_indices, 2:9)
    stopifnot(whole$requested_reads == 20, whole$rolling_reads == 8,
              whole$whole_saved_fraction == .6, whole$rolling_saved_fraction == .6)
    between <- (whole$whole_bytes + whole$rolling_bytes) / 2 / 1e9
    rolling <- planPgenGenotypeCache(reader, worklist, selected, cache_memory_gb = between)
    expect_identical_value(rolling$mode, "rolling")
    none <- planPgenGenotypeCache(reader, worklist, selected,
      cache_memory_gb = whole$rolling_bytes / 2 / 1e9)
    expect_identical_value(none$mode, "none")
    boundary <- planPgenGenotypeCache(reader, worklist, selected, mode = "whole",
      cache_memory_gb = whole$whole_bytes / 1e9)
    expect_identical_value(boundary$mode, "whole")
    sparse <- worklist[c(1L, 4L), ]
    sparse$window_start <- sparse$window_end <- c(100L, 200L)
    expect_identical_value(planPgenGenotypeCache(reader, sparse, selected)$mode, "none")
  })
  run_named_test("cache configuration rejects invalid modes and memory", {
    for (value in c(-1, NA_real_, Inf, NaN)) {
      expect_error_message(planPgenGenotypeCache(reader, worklist, selected,
        cache_memory_gb = value), "memory")
    }
    expect_error_message(planPgenGenotypeCache(reader, worklist, selected,
      job_memory_gb = 0), "memory")
    expect_error_message(planPgenGenotypeCache(reader, worklist, selected,
      mode = "bad"), "mode")
    expect_error_message(planPgenGenotypeCache(reader, worklist, selected,
      mode = "whole", cache_memory_gb = 1e-12), "budget")
    expect_identical_value(planPgenGenotypeCache(reader, worklist, selected,
      eligible = FALSE, reason = "additional genotypes")$mode, "none")
    expect_error_message(planPgenGenotypeCache(reader, worklist, selected,
      eligible = FALSE, mode = "whole", reason = "additional genotypes"), "additional genotypes")
  })
  run_named_test("cache plans report why sufficient statistics are disabled", {
    plan <- planPgenGenotypeCache(reader, worklist, selected, mode = "whole")
    stopifnot(is.list(plan$sufficient_statistics),
              !plan$sufficient_statistics$enabled,
              grepl("sample", plan$sufficient_statistics$reason))
  })
  reference <- function(indices, input_design = design) {
    mat <- calls[indices, selected, drop = FALSE]
    keep <- rowSums(!is.na(mat)) > 0L
    mat <- mat[keep, , drop = FALSE]
    indices <- indices[keep]
    means <- rowMeans(mat, na.rm = TRUE)
    missing <- which(is.na(mat), arr.ind = TRUE)
    if (nrow(missing)) mat[missing] <- means[missing[, 1L]]
    keep <- rowSums(round(mat, 0)) != 0
    mat <- mat[keep, , drop = FALSE]
    indices <- indices[keep]
    x <- t(mat - rowMeans(mat))
    x <- x - input_design %*% solve(crossprod(input_design)) %*% crossprod(input_design, x)
    dimnames(x) <- list(selected, reader$variants$variant_id[indices])
    list(x = x, indices = as.integer(indices))
  }
  run_named_test("whole cache prepares each physical variant once and matches R", {
    plan <- planPgenGenotypeCache(reader, worklist, selected, mode = "whole")
    cache <- createPgenGenotypeCache(reader, plan, selected, design)
    on.exit(clearPgenGenotypeCache(cache), add = TRUE)
    for (id in names(plan$window_indices)) {
      result <- getPgenGenotypeCacheWindow(cache, id)
      expected <- reference(plan$window_indices[[id]])
      stopifnot(isTRUE(all.equal(result$gt_hat, expected$x, tolerance = 1e-13)))
      expect_identical_value(result$variant_indices, expected$indices)
      stopifnot(!result$guard_required)
    }
    summary <- pgenGenotypeCacheSummary(cache)
    stopifnot(summary$variants_read == 8, summary$variants_prepared == 6,
              summary$excluded_variants == 2, summary$cache_hits > 0,
              summary$resident_bytes <= plan$budget_bytes)
    clearPgenGenotypeCache(cache)
    stopifnot(pgenGenotypeCacheSummary(cache)$resident_variants == 0)
    expect_error_message(getPgenGenotypeCacheWindow(cache, "A"), "cleared")
  })
  run_named_test("rolling cache reads only entering variants and rereads evicted variants", {
    plan <- planPgenGenotypeCache(reader, worklist, selected, mode = "rolling")
    cache <- createPgenGenotypeCache(reader, plan, selected, design)
    on.exit(clearPgenGenotypeCache(cache), add = TRUE)
    observed <- observe_reads({
      for (id in c("A", "B", "D", "A")) {
        result <- getPgenGenotypeCacheWindow(cache, id)
        expected <- reference(plan$window_indices[[id]])
        stopifnot(isTRUE(all.equal(result$gt_hat, expected$x, tolerance = 1e-13)))
        expect_identical_value(result$variant_indices, expected$indices)
      }
    })
    # The native retry removes all-missing physical row 5 only. Later reads
    # contain entering rows only; evicted ordinary rows 2 and 4 are read again.
    expect_identical_value(observed$indices,
      list(2:6, c(2L, 3L, 4L, 6L), 7L, 8:9, c(2L, 4L)), "native entering-variant reads")
    summary <- pgenGenotypeCacheSummary(cache)
    # A reads 5 rows (2 exclusions), B reads 7, D reads 8,9, then A rereads 2,4.
    stopifnot(summary$variants_read == 10, summary$cache_hits > 0,
              summary$variants_evicted > 0, summary$resident_variants == 3,
              summary$excluded_variants == 2)
  })
  run_named_test("cache flags covariate-explained columns without applying dense math", {
    explained_design <- cbind(variant = calls[2L, selected], intercept = 1)
    rownames(explained_design) <- selected
    plan <- planPgenGenotypeCache(reader, worklist, selected, mode = "whole")
    cache <- createPgenGenotypeCache(reader, plan, selected, explained_design)
    on.exit(clearPgenGenotypeCache(cache), add = TRUE)
    result <- getPgenGenotypeCacheWindow(cache, "A")
    stopifnot(result$guard_required,
              reader$variants$variant_id[[2L]] %in% result$covariate_explained_variants$variant_id)
    expected <- reference(plan$window_indices$A, explained_design)
    stopifnot(isTRUE(all.equal(result$gt_hat, expected$x, tolerance = 1e-13)))
  })
  run_named_test("indexed reader preserves physical order and returns empty imputed blocks", {
    result <- readPgenVariants(reader, c(7L, 2L, 5L, 4L), sample_ids = selected)
    expect_identical_value(result$variant_indices, c(7L, 2L, 4L))
    expect_identical_value(colnames(result$genotypes), selected)
    empty <- readPgenVariants(reader, 5L, sample_ids = selected)
    stopifnot(nrow(empty$genotypes) == 0, ncol(empty$genotypes) == length(selected),
              identical(empty$variant_indices, integer()))
  })
  run_named_test("cache context validation prevents sample or input reuse errors", {
    plan <- planPgenGenotypeCache(reader, worklist, selected, mode = "whole")
    expect_error_message(createPgenGenotypeCache(reader, plan, rev(selected), design), "samples")
    expect_error_message(createPgenGenotypeCache(reader, plan, selected,
      design[rev(seq_len(nrow(design))), , drop = FALSE]), "sample order")
    cache <- createPgenGenotypeCache(reader, plan, selected, design)
    on.exit(clearPgenGenotypeCache(cache), add = TRUE)
    expect_error_message(getPgenGenotypeCacheWindow(cache, "absent"), "absent")
    original_path <- reader$pgen_path
    reader$pgen_path <- paste0(original_path, ".changed")
    tryCatch({
      expect_error_message(getPgenGenotypeCacheWindow(cache, "A"), "identity")
      expect_error_message(createPgenGenotypeCache(reader, plan, selected, design), "PGEN")
    }, finally = reader$pgen_path <- original_path)
  })
  run_named_test("bounded blocks retain duplicate variant keys and physical window order", {
    large_files <- file.path(directory, c("large.pgen", "large.pvar", "large.psam"))
    base <- c(0, 1, 2, 0, 2, 1, 1, 0)
    large_calls <- t(vapply(seq_len(600L), function(i)
      base[((seq_len(8L) + i) %% 8L) + 1L], numeric(8L)))
    colnames(large_calls) <- samples
    large_calls[257L, ] <- NA_real_
    large_calls[511L, ] <- 0
    large_writer <- pgenlibr::NewPgenWriter(large_files[[1L]], 8L, 600L, nonref_flags = FALSE)
    tryCatch({
      for (i in seq_len(600L)) pgenlibr::AppendBiallelic(large_writer, large_calls[i, ])
    }, finally = pgenlibr::ClosePgenWriter(large_writer))
    large_positions <- seq_len(600L)
    large_positions[[300L]] <- 299L
    writeLines(c("#CHROM\tPOS\tID\tREF\tALT", paste("1", large_positions,
      paste0("v", seq_len(600L)), "A", "G", sep = "\t")), large_files[[2L]])
    writeLines(c("#IID", samples), large_files[[3L]])
    large_reader <- do.call(openPgenReader, as.list(large_files))
    on.exit(closePgenReader(large_reader), add = TRUE)
    large_windows <- data.frame(outcome_id = c("left", "right", "return"), chromosome = "1",
      window_start = c(1L, 201L, 1L), window_end = c(400L, 600L, 400L))
    plan <- planPgenGenotypeCache(large_reader, large_windows, selected, mode = "whole")
    cache <- createPgenGenotypeCache(large_reader, plan, selected, design)
    on.exit(clearPgenGenotypeCache(cache), add = TRUE)
    observed <- observe_reads(getPgenGenotypeCacheWindow(cache, "left"))
    result <- observed$value
    keep <- setdiff(seq_len(400L), 257L)
    mat <- large_calls[keep, selected, drop = FALSE]
    values <- t(mat - rowMeans(mat))
    expected <- values - design %*% solve(crossprod(design)) %*% crossprod(design, values)
    dimnames(expected) <- list(selected, large_reader$variants$variant_id[keep])
    stopifnot(isTRUE(all.equal(result$gt_hat, expected, tolerance = 1e-13)),
              max(lengths(observed$indices)) <= 256L,
              length(cache$blocks) == 3L,
              pgenGenotypeCacheSummary(cache)$variants_read == 600,
              pgenGenotypeCacheSummary(cache)$resident_variants == 598)
    expect_identical_value(result$variant_indices, as.integer(keep))
    stopifnot(sum(result$variant_names == "chr1_299_A_G") == 2L)
    # A disjoint/new portion was prepared with the union and needs no new reads.
    observed <- observe_reads(getPgenGenotypeCacheWindow(cache, "right"))
    stopifnot(length(observed$indices) == 0L)
  })
  ss_files <- file.path(directory, c("statistics.pgen", "statistics.pvar", "statistics.psam"))
  ss_samples <- paste0("n", seq_len(80L))
  ss_selected <- rev(ss_samples)
  set.seed(195L)
  ss_calls <- matrix(as.double(rbinom(40L * 80L, 2L, .3)), nrow = 40L,
    dimnames = list(NULL, ss_samples))
  ss_calls[2L, 7L] <- NA_real_
  ss_calls[12L, ] <- NA_real_
  ss_calls[18L, ] <- 0
  ss_writer <- pgenlibr::NewPgenWriter(ss_files[[1L]], 80L, 40L, nonref_flags = FALSE)
  tryCatch({
    for (i in seq_len(40L)) pgenlibr::AppendBiallelic(ss_writer, ss_calls[i, ])
  }, finally = pgenlibr::ClosePgenWriter(ss_writer))
  ss_positions <- seq_len(40L)
  ss_positions[[30L]] <- 29L
  writeLines(c("#CHROM\tPOS\tID\tREF\tALT", paste("1", ss_positions,
    paste0("v", seq_len(40L)), "A", "G", sep = "\t")), ss_files[[2L]])
  writeLines(c("#IID", ss_samples), ss_files[[3L]])
  ss_reader <- do.call(openPgenReader, as.list(ss_files))
  on.exit(closePgenReader(ss_reader), add = TRUE)
  ss_worklist <- data.frame(outcome_id = c("left", "right"), chromosome = "1",
    window_start = c(1L, 11L), window_end = c(30L, 40L))
  ss_design <- cbind(pc = rnorm(80L), intercept = 1)
  rownames(ss_design) <- ss_selected
  ss_plan <- planPgenGenotypeCache(ss_reader, ss_worklist, ss_selected, mode = "whole")
  run_named_test("sufficient statistics require overlap, sample count and peak memory", {
    stopifnot(ss_plan$sufficient_statistics$enabled,
              ss_plan$sufficient_statistics$estimated_bytes == 146240,
              ss_plan$sufficient_statistics$xtx_bytes == 12800,
              ss_plan$sufficient_statistics$working_bytes == 31040)
    # These limits include the resident cache, the shared Gram matrix, and
    # the largest native per-window matrix temporaries.
    boundary <- ss_plan
    boundary$budget_bytes <- 146240
    stopifnot(planPgenSufficientStatistics(boundary)$enabled)
    boundary$budget_bytes <- 146239
    stopifnot(!planPgenSufficientStatistics(boundary)$enabled)
    for (mode in c("rolling", "none")) {
      changed <- ss_plan
      changed$mode <- mode
      stopifnot(!planPgenSufficientStatistics(changed)$enabled)
    }
    changed <- ss_plan
    changed$sample_count <- 59L
    stopifnot(!planPgenSufficientStatistics(changed)$enabled)
    changed <- ss_plan
    changed$window_indices$right <- integer()
    stopifnot(!planPgenSufficientStatistics(changed)$enabled)
    changed <- ss_plan
    changed$max_window <- 10L
    changed$window_indices <- split(seq_len(40L), rep(seq_len(4L), each = 10L))
    stopifnot(!planPgenSufficientStatistics(changed)$enabled)
  })
  run_named_test("shared statistics match centered window crossproducts without extra reads", {
    cache <- createPgenGenotypeCache(ss_reader, ss_plan, ss_selected, ss_design)
    on.exit(clearPgenGenotypeCache(cache), add = TRUE)
    stopifnot(is.null(cache$sufficient_statistics))
    windows <- lapply(c("left", "right"), function(id) getPgenGenotypeCacheWindow(cache, id))
    observed <- observe_reads(lapply(windows, function(window)
      getPgenGenotypeCacheSufficientStatistics(cache, window)))
    stopifnot(length(observed$indices) == 0L)
    for (i in seq_along(windows)) {
      X <- windows[[i]]$gt_hat
      centered <- sweep(X, 2L, colMeans(X), "-")
      actual <- observed$value[[i]]
      stopifnot(isTRUE(all.equal(actual$XtX, crossprod(centered), tolerance = 1e-13)),
                isTRUE(all.equal(actual$X_colmeans, colMeans(X), tolerance = 1e-13)),
                identical(actual$n, 80L),
                identical(dimnames(actual$XtX), list(colnames(X), colnames(X))))
    }
    stopifnot(sum(colnames(windows[[1L]]$gt_hat) == "chr1_29_A_G") == 2L)
    summary <- pgenGenotypeCacheSummary(cache)
    stopifnot(summary$sufficient_statistics_enabled,
              summary$sufficient_statistics_builds == 1L,
              summary$sufficient_statistics_hits == 1L,
              summary$sufficient_statistics_bytes == 8 * 38^2,
              summary$variants_read == 40L)
    clearPgenGenotypeCache(cache)
    stopifnot(is.null(cache$sufficient_statistics),
              pgenGenotypeCacheSummary(cache)$sufficient_statistics_bytes == 0)
    expect_error_message(getPgenGenotypeCacheSufficientStatistics(cache, windows[[1L]]), "cleared")
  })
  run_named_test("shared statistics reject changed cache inputs and physical window indices", {
    cache <- createPgenGenotypeCache(ss_reader, ss_plan, ss_selected, ss_design)
    on.exit(clearPgenGenotypeCache(cache), add = TRUE)
    window <- getPgenGenotypeCacheWindow(cache, "left")
    changed <- window
    changed$variant_indices[[1L]] <- 40L
    expect_error_message(getPgenGenotypeCacheSufficientStatistics(cache, changed), "physical")
    original_path <- ss_reader$pgen_path
    ss_reader$pgen_path <- paste0(original_path, ".changed")
    tryCatch({
      expect_error_message(getPgenGenotypeCacheSufficientStatistics(cache, window), "identity")
    }, finally = ss_reader$pgen_path <- original_path)
    window$guard_required <- TRUE
    stopifnot(is.null(getPgenGenotypeCacheSufficientStatistics(cache, window)),
              is.null(cache$sufficient_statistics))
    cache$plan$sufficient_statistics$enabled <- FALSE
    window$guard_required <- FALSE
    stopifnot(is.null(getPgenGenotypeCacheSufficientStatistics(cache, window)),
              is.null(cache$sufficient_statistics))
  })
  message("All prepared PGEN cache tests passed.")
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  benchmark_library <- Sys.getenv("BENCH_LIBRARY", "")
  if (nzchar(benchmark_library)) .libPaths(c(benchmark_library, .libPaths()))
  args <- commandArgs(trailingOnly = TRUE)
  repo_dir <- normalizePath(if (length(args)) args[[1L]] else getwd(), mustWork = TRUE)
  run_pgen_cache_tests(repo_dir)
}
