# Prepared cache state tests do not require an installed PGEN reader.
run_pgen_cache_state_tests <- function(repo_dir) {
  source(file.path(repo_dir, "tests", "test_helpers.R"), local = TRUE)
  functions_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH", file.path(repo_dir, "R", "utils"))
  source(file.path(functions_path, "PgenCacheFunctions.R"), local = TRUE)
  samples <- c("s3", "s1")
  reader <- new.env(parent = emptyenv())
  class(reader) <- "PgenReader"
  reader$closed <- FALSE
  reader$sample_ids <- samples
  reader$pgen_path <- "prepared-cache-state.pgen"
  reader$variants <- data.frame(variant_id = c("unused", "duplicate", "v3", "unused",
                                              "duplicate", "unused", "unused", "v8"))
  cache <- new.env(parent = emptyenv())
  class(cache) <- "PgenGenotypeCache"
  cache$reader <- reader
  cache$pgen_path <- reader$pgen_path
  cache$sample_ids <- samples
  cache$mode <- "whole"
  cache$initialized <- TRUE
  cache$cleared <- FALSE
  cache$cache_hits <- cache$exclusion_hits <- cache$variants_evicted <- 0L
  cache$excluded_indices <- integer()
  cache$plan <- list(window_indices = list(left = c(2L, 3L, 5L), right = c(3L, 5L, 8L),
                                           empty = integer()))
  cache$blocks <- list(
    list(indices = c(5L, 2L), gt_hat = matrix(c(5, -5, 2, -2), nrow = 2L),
         cohort_maf = c(.375, .125), cohort_an = c(6, 8), cohort_mac = c(2, 1),
         flagged_indices = integer(),
         diagnostics = .pgenCacheEmptyDiagnostics()),
    list(indices = c(8L, 3L), gt_hat = matrix(c(8, -8, 3, -3), nrow = 2L),
         cohort_maf = c(.5, .25), cohort_an = c(4, 6), cohort_mac = c(2, 1.5),
         flagged_indices = integer(),
         diagnostics = .pgenCacheEmptyDiagnostics()))
  run_named_test("cached MAF follows physical predictors with repeated labels", {
    left <- getPgenGenotypeCacheWindow(cache, "left")
    expect_identical_value(left$cohort_maf, c(.125, .25, .375), "left physical MAF order")
    expect_identical_value(left$variant_names, c("duplicate", "v3", "duplicate"))
    expect_identical_value(left$cohort_an, c(8, 6, 6), "left observed AN order")
    expect_identical_value(left$cohort_mac, c(1, 1.5, 2), "left observed MAC order")
    right <- getPgenGenotypeCacheWindow(cache, "right")
    expect_identical_value(right$cohort_maf, c(.25, .375, .5), "right physical MAF order")
    expect_identical_value(right$cohort_an, c(6, 6, 4), "right observed AN order")
    expect_identical_value(right$cohort_mac, c(1.5, 2, 2), "right observed MAC order")
  })
  run_named_test("rolling eviction keeps MAF aligned with resident genotype columns", {
    cache$mode <- "rolling"
    left <- getPgenGenotypeCacheWindow(cache, "left")
    expect_identical_value(left$cohort_maf, c(.125, .25, .375), "evicted window MAF")
    expect_identical_value(cache$blocks[[2L]]$cohort_maf, .25, "retained block MAF")
    expect_identical_value(cache$blocks[[2L]]$cohort_an, 6, "retained block AN")
    expect_identical_value(cache$blocks[[2L]]$cohort_mac, 1.5, "retained block MAC")
    empty <- getPgenGenotypeCacheWindow(cache, "empty")
    expect_identical_value(empty$cohort_maf, numeric(), "empty window MAF")
    expect_identical_value(empty$cohort_an, numeric(), "empty window AN")
    expect_identical_value(empty$cohort_mac, numeric(), "empty window MAC")
  })
  invisible(TRUE)
}

# Real PGEN files test cache selection, physical variant mapping, and reuse.
run_pgen_cache_tests <- function(repo_dir) {
  run_pgen_cache_state_tests(repo_dir)
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
  calls[7L, ] <- c(2, 0, 2, 1, 2, 0, 2, 1)
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
    observer$raw_indices <- list()
    tracer <- substitute({
      .cache_read_observer <- captured
      if (meanimpute) {
        .cache_read_observer$indices[[length(.cache_read_observer$indices) + 1L]] <- variant_subset
      } else {
        .cache_read_observer$raw_indices[[length(.cache_read_observer$raw_indices) + 1L]] <- variant_subset
      }
    }, list(captured = observer))
    suppressMessages(trace("ReadList", tracer = tracer,
      where = asNamespace("pgenlibr"), print = FALSE))
    on.exit(suppressMessages(untrace("ReadList", where = asNamespace("pgenlibr"))), add = TRUE)
    value <- force(code)
    list(value = value, indices = observer$indices, raw_indices = observer$raw_indices)
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
    alt_af <- rowMeans(mat) / 2
    observed <- calls[indices, selected, drop = FALSE]
    an <- 2 * rowSums(!is.na(observed))
    ac <- rowSums(observed, na.rm = TRUE)
    list(x = x, indices = as.integer(indices), maf = unname(pmin(alt_af, 1 - alt_af)),
         an = unname(as.double(an)), mac = unname(as.double(pmin(ac, an - ac))))
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
      stopifnot(isTRUE(all.equal(result$cohort_maf, expected$maf, tolerance = 1e-15)))
      expect_identical_value(result$cohort_an, expected$an, "whole-cache observed AN")
      expect_identical_value(result$cohort_mac, expected$mac, "whole-cache observed MAC")
      stopifnot(!result$guard_required)
    }
    first <- getPgenGenotypeCacheWindow(cache, "A")
    stopifnot(isTRUE(all.equal(first$cohort_maf[first$variant_indices == 2L], 5 / 12,
                              tolerance = 1e-15)))
    high_alt <- getPgenGenotypeCacheWindow(cache, "B")
    stopifnot(isTRUE(all.equal(high_alt$cohort_maf[high_alt$variant_indices == 7L], 1 / 6,
                              tolerance = 1e-15)))
    expect_identical_value(high_alt$cohort_an[high_alt$variant_indices == 7L], 12,
                           "ALT-major AN")
    expect_identical_value(high_alt$cohort_mac[high_alt$variant_indices == 7L], 2,
                           "ALT-major minor REF count")
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
        stopifnot(isTRUE(all.equal(result$cohort_maf, expected$maf, tolerance = 1e-15)))
        expect_identical_value(result$cohort_an, expected$an, "rolling-cache observed AN")
        expect_identical_value(result$cohort_mac, expected$mac, "rolling-cache observed MAC")
      }
    })
    # The native retry removes all-missing physical row 5 only. Later reads
    # contain entering rows only; evicted ordinary rows 2 and 4 are read again.
    expect_identical_value(observed$indices,
      list(2:6, c(2L, 3L, 4L, 6L), 7L, 8:9, c(2L, 4L)), "native entering-variant reads")
    expect_identical_value(observed$raw_indices,
      list(2:6, 7L, 8:9, c(2L, 4L)), "observed counts read only entering variants")
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
              max(lengths(observed$raw_indices)) <= 256L,
              length(cache$blocks) == 3L,
              pgenGenotypeCacheSummary(cache)$variants_read == 600,
              pgenGenotypeCacheSummary(cache)$resident_variants == 598)
    expect_identical_value(result$variant_indices, as.integer(keep))
    stopifnot(sum(result$variant_names == "chr1_299_A_G") == 2L)
    stopifnot(isTRUE(all.equal(result$cohort_maf[result$variant_names == "chr1_299_A_G"],
                              c(5 / 12, 1 / 3), tolerance = 1e-15)))
    # A disjoint/new portion was prepared with the union and needs no new reads.
    observed <- observe_reads(getPgenGenotypeCacheWindow(cache, "right"))
    stopifnot(length(observed$indices) == 0L)
    stopifnot(length(observed$raw_indices) == 0L)
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
