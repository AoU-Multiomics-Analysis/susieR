# These checks use real PGEN calls. A change to the imputation cohort, variant
# mapping, or dosage handling must change the returned values and fail here.
run_pgen_native_imputation_tests <- function(repo_dir) {
  source(file.path(repo_dir, "tests", "test_helpers.R"), local = TRUE)
  functions_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH",
                               unset = file.path(repo_dir, "R", "utils"))
  for (filename in c("SusieFunctions.R", "PgenFunctions.R", "PgenGeneFunctions.R")) {
    source(file.path(functions_path, filename), local = TRUE)
  }
  directory <- tempfile("pgen-native-imputation-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  samples <- c("s6", "001", "s2", "s3", "NA", "s1")
  ids <- paste0("chr1_", 101:106, "_A_G")
  calls <- rbind(
    c(2, 0, NA, 1, 2, 2),
    c(0, 2, 1, NA, 2, 2),
    c(0, 1, 2, 0, 1, 2),
    rep(NA_real_, 6L),
    c(1, NA, 2, 2, 0, 1),
    c(0, NA, 0, 0, 0, 0)
  )
  dimnames(calls) <- list(ids, samples)
  files <- file.path(directory, c("input.pgen", "input.pvar", "input.psam"))
  writer <- pgenlibr::NewPgenWriter(files[[1L]], length(samples), nrow(calls),
                                   nonref_flags = FALSE)
  tryCatch({
    for (i in seq_len(nrow(calls))) pgenlibr::AppendBiallelic(writer, calls[i, ])
  }, finally = pgenlibr::ClosePgenWriter(writer))
  writeLines(c("#CHROM\tPOS\tID\tREF\tALT",
    paste("1", 101:106, ids, "A", "G", sep = "\t")), files[[2L]])
  writeLines(c("#IID", samples), files[[3L]])
  reader <- do.call(openPgenReader, as.list(files))
  on.exit(closePgenReader(reader), add = TRUE)

  # Keep the native function real. The observer records its inputs and can
  # inject an error at its boundary to check temporary-handle cleanup.
  observe_native <- function(code, fail = FALSE) {
    observer <- new.env(parent = emptyenv())
    observer$calls <- list()
    observer$hardcall_calls <- 0L
    observer$fail <- fail
    tracer <- substitute({
      .native_imputation_observer <- .observer
      .native_imputation_observer$calls[[length(.native_imputation_observer$calls) + 1L]] <- list(
        variants = variant_subset, meanimpute = meanimpute, handle = pgen)
      if (.native_imputation_observer$fail) {
        stop("injected native imputation error", call. = FALSE)
      }
    }, list(.observer = observer))
    suppressMessages(trace("ReadList", tracer = tracer,
      where = asNamespace("pgenlibr"), print = FALSE))
    hardcall_tracer <- substitute({
      .native_imputation_observer <- .observer
      .native_imputation_observer$hardcall_calls <-
        .native_imputation_observer$hardcall_calls + 1L
    }, list(.observer = observer))
    suppressMessages(trace("ReadIntList", tracer = hardcall_tracer,
      where = asNamespace("pgenlibr"), print = FALSE))
    on.exit({
      suppressMessages(untrace("ReadList", where = asNamespace("pgenlibr")))
      suppressMessages(untrace("ReadIntList", where = asNamespace("pgenlibr")))
    }, add = TRUE)
    value <- tryCatch(force(code), error = identity)
    list(value = value, calls = observer$calls, hardcall_calls = observer$hardcall_calls)
  }
  selected <- c("s3", "s2", "001")
  rows <- c(2L, 3L, 1L, 5L)
  input <- calls[rows, selected, drop = FALSE]
  expected <- rbind(c(1.5, 1, 2), c(0, 2, 1), c(1, 0.5, 0), c(2, 2, 2))
  dimnames(expected) <- dimnames(input)
  run_named_test("native imputation uses exact selected samples and preserves requested order", {
    observed <- observe_native(imputePgenHardcalls(reader, input))
    if (inherits(observed$value, "error")) stop(observed$value)
    expect_identical_value(observed$value, expected, "selected-sample means")
    stopifnot(length(observed$calls) > 0L,
      all(vapply(observed$calls, function(call) isTRUE(call$meanimpute), logical(1))))
    expect_identical_value(sort(unique(unlist(lapply(observed$calls, `[[`, "variants")))),
                           c(1L, 2L, 5L), "only variants with missing calls are reread")
    for (call in observed$calls) {
      expect_error_message(pgenlibr::GetVariantCt(call$handle), "closed")
    }
    expect_identical_value(readPgenWindow(reader, "1", 101L, 106L, as_double = TRUE),
                           calls, "full-cohort reader after imputation")
  })
  run_named_test("native imputation handles the minus-one missing policy", {
    sentinel <- input
    sentinel[is.na(sentinel)] <- -1
    expect_identical_value(imputePgenHardcalls(reader, sentinel), expected,
                           "minus-one imputation")
  })
  run_named_test("native imputation retains relative inputs after a directory change", {
    original_directory <- getwd()
    setwd(directory)
    relative_reader <- openPgenReader("input.pgen", "input.pvar", "input.psam")
    tryCatch({
      setwd(dirname(directory))
      expect_identical_value(readPgenWindow(relative_reader, "1", 101L, 106L,
                                            as_double = TRUE), calls)
      expect_identical_value(imputePgenHardcalls(relative_reader, input), expected,
                             "imputation after changing the working directory")
    }, finally = {
      setwd(original_directory)
      closePgenReader(relative_reader)
    })
  })
  run_named_test("native imputation preserves one-row and one-column matrix shapes", {
    expect_identical_value(imputePgenHardcalls(reader, input[1L, , drop = FALSE]),
                           expected[1L, , drop = FALSE], "one variant")
    one_sample <- calls[c(1L, 3L), "001", drop = FALSE]
    expect_identical_value(imputePgenHardcalls(reader, one_sample), one_sample,
                           "one sample without missing values")
  })
  run_named_test("complete genotypes skip native imputation", {
    complete <- calls[3L, selected, drop = FALSE]
    observed <- observe_native(imputePgenHardcalls(reader, complete))
    if (inherits(observed$value, "error")) stop(observed$value)
    expect_identical_value(observed$value, complete, "complete genotypes")
    stopifnot(length(observed$calls) == 0L)
  })
  run_named_test("all-missing variants are removed before native imputation", {
    for (policy in c("NA", "-1")) {
      missing <- calls[c(4L, 2L, 3L), selected, drop = FALSE]
      if (policy == "-1") missing[is.na(missing)] <- -1
      expected_missing <- rbind(c(1.5, 1, 2), c(0, 2, 1))
      dimnames(expected_missing) <- list(ids[c(2L, 3L)], selected)
      observed <- observe_native(imputePgenHardcalls(reader, missing))
      if (inherits(observed$value, "error")) stop(observed$value)
      expect_identical_value(observed$value, expected_missing,
                             paste(policy, "all-missing variant removed"))
      expect_identical_value(unique(unlist(lapply(observed$calls, `[[`, "variants"))),
                             2L, "native reader receives only nonempty missing rows")
      all_missing <- missing[1L, , drop = FALSE]
      empty <- matrix(numeric(), nrow = 0L, ncol = length(selected),
                      dimnames = list(NULL, selected))
      observed <- observe_native(imputePgenHardcalls(reader, all_missing))
      if (inherits(observed$value, "error")) stop(observed$value)
      expect_identical_value(observed$value, empty, "no usable variants remain")
      stopifnot(length(observed$calls) == 0L)
    }
  })
  run_named_test("integer-valued stored dosages match the hardcall calculation", {
    dosage_files <- file.path(directory, c("dosage.pgen", "dosage.pvar", "dosage.psam"))
    dosage_writer <- pgenlibr::NewPgenWriter(dosage_files[[1L]], 6L, 1L,
      nonref_flags = FALSE, dosage_present = TRUE)
    tryCatch(pgenlibr::AppendDosages(dosage_writer, c(0, 1, 2, 1, NA, 2)),
             finally = pgenlibr::ClosePgenWriter(dosage_writer))
    writeLines(c("#CHROM\tPOS\tID\tREF\tALT", "1\t201\tdosage\tA\tG"),
               dosage_files[[2L]])
    writeLines(c("#IID", samples), dosage_files[[3L]])
    dosage_reader <- do.call(openPgenReader, as.list(dosage_files))
    tryCatch({
      hardcalls <- readPgenWindow(dosage_reader, "1", 201L, 201L, as_double = TRUE)
      hardcalls <- hardcalls[, c("NA", "s6", "s2", "s3"), drop = FALSE]
      expected_dosage <- matrix(c(1, 0, 2, 1), nrow = 1L,
                                dimnames = dimnames(hardcalls))
      expect_identical_value(imputePgenHardcalls(dosage_reader, hardcalls),
                             expected_dosage, "integer-valued dosage means")
    }, finally = closePgenReader(dosage_reader))
  })
  run_named_test("multiple native read blocks preserve filtered variant and sample order", {
    large_files <- file.path(directory, c("large.pgen", "large.pvar", "large.psam"))
    large_samples <- c("s7", "001", "s6", "s3", "NA", "s5", "s4")
    large_selected <- large_samples[c(7L, 2L, 5L, 1L, 4L)]
    patterns <- rbind(c(0, 1, 2, NA, 2, 0, 1),
                      c(2, NA, 0, 1, 0, 2, 1),
                      c(NA, 0, 2, 2, 1, 0, 0))
    pattern_index <- rep(1:3, length.out = 306L)
    large_calls <- patterns[pattern_index, , drop = FALSE]
    large_ids <- paste0("chr1_", 1001:1306, "_A_G")
    dimnames(large_calls) <- list(large_ids, large_samples)
    # The selected cohort has no calls in these rows even though other PGEN
    # samples do. Reversing the rows places them before and across a read block.
    all_missing_rows <- c(306L, 50L, 49L)
    large_calls[all_missing_rows, large_selected] <- NA_real_
    large_writer <- pgenlibr::NewPgenWriter(large_files[[1L]], 7L, 306L,
                                           nonref_flags = FALSE)
    tryCatch({
      for (i in seq_len(nrow(large_calls))) {
        pgenlibr::AppendBiallelic(large_writer, large_calls[i, ])
      }
    }, finally = pgenlibr::ClosePgenWriter(large_writer))
    writeLines(c("#CHROM\tPOS\tID\tREF\tALT",
      paste("1", 1001:1306, large_ids, "A", "G", sep = "\t")), large_files[[2L]])
    writeLines(c("#IID", large_samples), large_files[[3L]])
    large_reader <- do.call(openPgenReader, as.list(large_files))
    tryCatch({
      feature_order <- rev(seq_len(nrow(large_calls)))
      input_large <- readPgenWindow(large_reader, "1", 1001L, 1306L, as_double = TRUE)
      input_large <- input_large[feature_order, large_selected, drop = FALSE]
      retained_order <- feature_order[!feature_order %in% all_missing_rows]
      # Literal means over the five selected samples: 1, 1, and 3/4.
      expected_patterns <- rbind(c(1, 1, 2, 0, 1),
                                 c(1, 1, 0, 2, 1),
                                 c(0, 0, 1, 0.75, 2))
      expected_large <- expected_patterns[pattern_index[retained_order], , drop = FALSE]
      dimnames(expected_large) <- list(large_ids[retained_order], large_selected)
      observed <- observe_native(imputePgenHardcalls(large_reader, input_large))
      if (inherits(observed$value, "error")) stop(observed$value)
      expect_identical_value(observed$value, expected_large, "values across native read blocks")
      stopifnot(length(observed$calls) > 1L)
      expect_identical_value(unname(unlist(lapply(observed$calls, `[[`, "variants"))),
                             as.integer(retained_order), "all retained native read indices")
    }, finally = closePgenReader(large_reader))
  })
  run_named_test("ambiguous variant names fail before reading the wrong genotype row", {
    duplicate_files <- file.path(directory, c("duplicate.pgen", "duplicate.pvar", "duplicate.psam"))
    duplicate_writer <- pgenlibr::NewPgenWriter(duplicate_files[[1L]], 6L, 2L,
                                              nonref_flags = FALSE)
    tryCatch({
      pgenlibr::AppendBiallelic(duplicate_writer, calls[1L, ])
      pgenlibr::AppendBiallelic(duplicate_writer, calls[2L, ])
    }, finally = pgenlibr::ClosePgenWriter(duplicate_writer))
    writeLines(c("#CHROM\tPOS\tID\tREF\tALT", "1\t301\ta\tA\tG", "1\t301\tb\tA\tG"),
               duplicate_files[[2L]])
    writeLines(c("#IID", samples), duplicate_files[[3L]])
    duplicate_reader <- do.call(openPgenReader, as.list(duplicate_files))
    tryCatch({
      duplicated <- readPgenWindow(duplicate_reader, "1", 301L, 301L, as_double = TRUE)
      duplicated <- duplicated[, selected, drop = FALSE]
      expect_error_message(imputePgenHardcalls(duplicate_reader, duplicated), "ambiguous")
    }, finally = closePgenReader(duplicate_reader))
  })
  run_named_test("native imputation closes temporary handles after an error", {
    observed <- observe_native(imputePgenHardcalls(reader, input), fail = TRUE)
    stopifnot(inherits(observed$value, "error"), length(observed$calls) == 1L)
    expect_error_message(stop(observed$value), "injected native imputation error")
    expect_error_message(pgenlibr::GetVariantCt(observed$calls[[1L]]$handle), "closed")
    expect_identical_value(readPgenWindow(reader, "1", 101L, 106L, as_double = TRUE),
                           calls, "source reader remains open after an error")
  })
  run_named_test("first native read imputes directly without a preceding hardcall read", {
    expected_first <- rbind(c(1, 0.5, 0), c(1.5, 1, 2), c(0, 2, 1))
    dimnames(expected_first) <- list(ids[1:3], selected)
    observed <- observe_native(readPgenWindow(reader, "1", 101L, 103L,
      meanimpute = TRUE, sample_ids = selected))
    if (inherits(observed$value, "error")) stop(observed$value)
    expect_identical_value(observed$value, expected_first, "direct native imputation")
    stopifnot(typeof(observed$value) == "double", observed$hardcall_calls == 0L,
      length(observed$calls) == 1L, isTRUE(observed$calls[[1L]]$meanimpute))
  })
  run_named_test("first native read is independent of the missing-code output policy", {
    expected_first <- rbind(c(1, 0.5, 0), c(1.5, 1, 2), c(0, 2, 1))
    dimnames(expected_first) <- list(ids[1:3], selected)
    for (policy in c("NA", "-1")) {
      actual <- readPgenWindow(reader, "1", 101L, 103L,
        as_double = FALSE, missing_code = policy, meanimpute = TRUE,
        sample_ids = selected)
      expect_identical_value(actual, expected_first, paste(policy, "native reader output"))
    }
  })
  run_named_test("first native read removes entirely missing variants and retries", {
    expected_first <- rbind(c(1, 0.5, 0), c(1.5, 1, 2), c(0, 2, 1),
                            c(2, 2, 2), c(0, 0, 0))
    dimnames(expected_first) <- list(ids[c(1L, 2L, 3L, 5L, 6L)], selected)
    observed <- observe_native(readPgenWindow(reader, "1", 101L, 106L,
      meanimpute = TRUE, sample_ids = selected))
    if (inherits(observed$value, "error")) stop(observed$value)
    expect_identical_value(observed$value, expected_first, "all-missing variant omitted")
    stopifnot(observed$hardcall_calls == 0L, length(observed$calls) == 2L)
    expect_identical_value(observed$calls[[2L]]$variants,
                           c(1L, 2L, 3L, 5L, 6L), "retry retains the other global indices")
    expect_error_message(readPgenWindow(reader, "1", 104L, 104L,
      meanimpute = TRUE, sample_ids = selected), "No variant with alt alleles")
    expect_identical_value(readPgenWindow(reader, "1", 101L, 106L, as_double = TRUE),
                           calls, "original reader after all-missing handling")
  })
  run_named_test("first native read propagates errors and closes its temporary handle", {
    observed <- observe_native(readPgenWindow(reader, "1", 101L, 103L,
      meanimpute = TRUE, sample_ids = selected), fail = TRUE)
    stopifnot(inherits(observed$value, "error"), length(observed$calls) == 1L,
              observed$hardcall_calls == 0L)
    expect_error_message(stop(observed$value), "injected native imputation error")
    expect_error_message(pgenlibr::GetVariantCt(observed$calls[[1L]]$handle), "closed")
    expect_identical_value(readPgenWindow(reader, "1", 101L, 106L, as_double = TRUE),
                           calls, "original reader remains open after a direct native error")
  })
  run_named_test("first native read preserves duplicate variant keys in PGEN file order", {
    duplicate_files <- file.path(directory, c("native-duplicate.pgen",
      "native-duplicate.pvar", "native-duplicate.psam"))
    duplicate_writer <- pgenlibr::NewPgenWriter(duplicate_files[[1L]], 6L, 3L,
                                              nonref_flags = FALSE)
    tryCatch({
      for (i in c(1L, 3L, 2L)) pgenlibr::AppendBiallelic(duplicate_writer, calls[i, ])
    }, finally = pgenlibr::ClosePgenWriter(duplicate_writer))
    writeLines(c("#CHROM\tPOS\tID\tREF\tALT", "1\t301\ta\tA\tG",
                 "1\t300\tmiddle\tA\tG", "1\t301\tb\tA\tG"), duplicate_files[[2L]])
    writeLines(c("#IID", samples), duplicate_files[[3L]])
    duplicate_reader <- do.call(openPgenReader, as.list(duplicate_files))
    tryCatch({
      expected_duplicate <- rbind(c(1, 0.5, 0), c(0, 2, 1), c(1.5, 1, 2))
      dimnames(expected_duplicate) <- list(c("chr1_301_A_G", "chr1_300_A_G",
                                            "chr1_301_A_G"), selected)
      observed <- observe_native(readPgenWindow(duplicate_reader, "1", 300L, 301L,
        meanimpute = TRUE, sample_ids = selected))
      if (inherits(observed$value, "error")) stop(observed$value)
      expect_identical_value(observed$value, expected_duplicate, "duplicate-key file order")
      stopifnot(observed$hardcall_calls == 0L, length(observed$calls) == 1L)
      expect_identical_value(readPgenWindow(duplicate_reader, "1", 301L, 301L,
        meanimpute = TRUE, sample_ids = selected), expected_duplicate[c(1L, 3L), , drop = FALSE],
        "duplicate-position window")
    }, finally = closePgenReader(duplicate_reader))
  })
  run_named_test("first native read validates sample subsets against the open reader", {
    subset_reader <- openPgenReader(files[[1L]], files[[2L]], files[[3L]], sample_ids = selected)
    tryCatch({
      requested <- rev(selected)
      expected_subset <- rbind(c(0, 0.5, 1), c(2, 1, 1.5), c(1, 2, 0))
      dimnames(expected_subset) <- list(ids[1:3], requested)
      expect_identical_value(readPgenWindow(subset_reader, "1", 101L, 103L,
        meanimpute = TRUE, sample_ids = requested), expected_subset, "reader sample subset")
      expect_error_message(readPgenWindow(subset_reader, "1", 101L, 103L,
        meanimpute = TRUE, sample_ids = "s6"), "sample_ids")
      for (bad in list(character(), NA_character_, c("001", "001"), 1)) {
        expect_error_message(readPgenWindow(subset_reader, "1", 101L, 103L,
          meanimpute = TRUE, sample_ids = bad), "sample_ids")
      }
    }, finally = closePgenReader(subset_reader))
  })
  message("All native PGEN imputation tests passed.")
  invisible(TRUE)
}
