# Small real PGEN fixtures test the public reader contract.
run_pgen_reader_tests <- function(repo_dir) {
  source(file.path(repo_dir, "tests", "test_helpers.R"), local = TRUE)
  functions_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH",
                               unset = file.path(repo_dir, "R", "utils"))
  source(file.path(functions_path, "PgenFunctions.R"), local = TRUE)
  if (!exists("openPgenReader", mode = "function")) {
    stop("Missing reader function: openPgenReader", call. = FALSE)
  }
  fixture_dir <- tempfile("pgen-reader-tests-")
  dir.create(fixture_dir)
  on.exit(unlink(fixture_dir, recursive = TRUE), add = TRUE)
  paths <- file.path(fixture_dir, c("small.pgen", "small.pvar", "small.psam"))
  sample_ids <- c("s3", "001", "NA", "s2", "s1", "s4")
  hardcalls <- rbind(
    c(2L, 0L, 2L, 1L, 0L, 1L),
    c(1L, 1L, 0L, 2L, 0L, 0L),
    c(0L, 1L, 2L, NA_integer_, 1L, 0L),
    c(0L, 0L, 0L, 0L, 0L, 0L),
    c(0L, NA_integer_, 0L, 0L, 0L, 0L),
    c(1L, NA_integer_, 0L, 0L, 0L, 0L),
    c(2L, NA_integer_, 0L, 0L, 0L, 0L),
    rep(NA_integer_, 6L),
    c(0L, 1L, 2L, 0L, 1L, 2L),
    c(2L, 0L, 1L, 2L, 0L, 1L),
    c(2L, 2L, 1L, 0L, NA_integer_, 1L),
    c(0L, 2L, 1L, 1L, 0L, 2L)
  )
  writer <- pgenlibr::NewPgenWriter(paths[[1L]], 6L, 12L, nonref_flags = FALSE)
  tryCatch({
    for (i in seq_len(nrow(hardcalls))) {
      pgenlibr::AppendBiallelic(writer, hardcalls[i, ])
    }
  }, finally = pgenlibr::ClosePgenWriter(writer))
  # Unsorted positions and duplicate positions must retain their file order.
  pvar_lines <- c(
    "##fileformat=PVARv1.2", "##fixture=reader",
    "#CHROM\tPOS\tID\tREF\tALT",
    "2\t90\tother\tA\tG", "1\t90\tbefore\tA\tG",
    "chr1\t100\tleft\tA\tG", "1\t110\tduplicate\tA\tT",
    "CHR1\t102\tzero_missing\tA\tG", "1\t103\tone_missing\tA\tG",
    "chr1\t104\ttwo_missing\tA\tG", "1\t105\tall_missing\tA\tG",
    "chr1\t106\tmixed\tA\tG", "1\t107\tNA\tNA\tG",
    "1\t110\tright\tA\tG", "1\t111\tafter\tA\tG"
  )
  psam_lines <- c("#FID\tIID\tSEX", paste("0", sample_ids, "0", sep = "\t"))
  writeLines(pvar_lines, paths[[2L]])
  writeLines(psam_lines, paths[[3L]])
  open_reader <- function(ids = NULL) {
    openPgenReader(paths[[1L]], paths[[2L]], paths[[3L]], sample_ids = ids)
  }
  reader <- open_reader()
  on.exit(closePgenReader(reader), add = TRUE)
  expected <- hardcalls[3:11, , drop = FALSE]
  dimnames(expected) <- list(c(
    "chr1_100_A_G", "chr1_110_A_T", "chr1_102_A_G", "chr1_103_A_G",
    "chr1_104_A_G", "chr1_105_A_G", "chr1_106_A_G", "chr1_107_NA_G",
    "chr1_110_A_G"
  ), sample_ids)
  run_named_test("global indices, exact integer hardcalls and missingness", {
    stopifnot(is.environment(reader), inherits(reader, "PgenReader"),
              tibble::is_tibble(reader$variants))
    expect_identical_value(reader$variants$pgen_index, 1:12)
    expect_identical_value(reader$variants$chromosome, c("2", rep("1", 11L)))
    expect_identical_value(reader$variants$original_id[[10L]], "NA")
    expect_identical_value(reader$variants$REF[[10L]], "NA")
    expect_identical_value(reader$variants$ALT[[10L]], "G")
    expect_identical_value(reader$variants$variant_id[3:11], rownames(expected))
    stopifnot(is.numeric(reader$variants$position))
    expect_identical_value(reader$sample_ids, sample_ids)
    expect_identical_value(reader$raw_sample_ids, sample_ids)
    expect_identical_value(reader$raw_sample_count, 6L)
    expect_identical_value(readPgenWindow(reader, "1", 100, 110), expected)
  })
  run_named_test("inclusive bounds, duplicate positions and file order", {
    expect_identical_value(readPgenWindow(reader, "chr1", 110, 110),
                           expected[c(2L, 9L), , drop = FALSE])
    expect_identical_value(readPgenWindow(reader, "CHRchr1", 100, 100),
                           expected[1L, , drop = FALSE])
    other <- hardcalls[1L, , drop = FALSE]
    dimnames(other) <- list("chr2_90_A_G", sample_ids)
    expect_identical_value(readPgenWindow(reader, 2, 90, 90), other)
  })
  run_named_test("double storage and explicit missing-code policies", {
    double <- expected
    storage.mode(double) <- "double"
    expect_identical_value(readPgenWindow(reader, "1", 100, 110, TRUE), double)
    legacy <- expected
    legacy[is.na(legacy)] <- -1L
    expect_identical_value(readPgenWindow(reader, "1", 100, 110,
                                         missing_code = "-1"), legacy)
    storage.mode(legacy) <- "double"
    expect_identical_value(readPgenWindow(reader, "1", 100, 110, TRUE, "-1"), legacy)
    expect_error_message(readPgenWindow(reader, "1", 100, 110,
                                        missing_code = "invalid"), "missing_code")
    expect_error_message(readPgenWindow(reader, "1", 100, 110,
                                        missing_code = NA_character_), "missing_code")
    expect_error_message(readPgenWindow(reader, "1", 100, 110,
                                        missing_code = c("NA", "-1")), "missing_code")
    expect_error_message(readPgenWindow(reader, "1", 100, 110,
                                        missing_code = -1), "missing_code")
    expect_error_message(readPgenWindow(reader, "1", 100, 110, NA), "as_double")
    expect_error_message(readPgenWindow(reader, "1", 100, 110, 1), "as_double")
    expect_error_message(readPgenWindow(reader, "1", 100, 110,
                                        c(TRUE, FALSE)), "as_double")
  })
  run_named_test("sample subsets, literal IDs and single dimensions", {
    subset_ids <- c("s1", "001", "NA")
    subset_reader <- open_reader(subset_ids)
    tryCatch({
      expect_identical_value(subset_reader$sample_ids, subset_ids)
      expect_identical_value(subset_reader$raw_sample_ids, sample_ids)
      expect_identical_value(subset_reader$raw_sample_count, 6L)
      expect_identical_value(readPgenWindow(subset_reader, "1", 100, 110),
                             expected[, subset_ids, drop = FALSE])
    }, finally = closePgenReader(subset_reader))
    single <- open_reader("001")
    tryCatch({
      expect_identical_value(readPgenWindow(single, "1", 100, 100),
                             expected[1L, "001", drop = FALSE])
      expect_identical_value(readPgenWindow(single, "1", 100, 110),
                             expected[, "001", drop = FALSE])
    }, finally = closePgenReader(single))
    expect_error_message(open_reader("absent"), "absent")
    expect_error_message(open_reader(c("s1", "s1")), "duplicate")
    expect_error_message(open_reader(character()), "sample_ids")
    expect_error_message(open_reader(NA_character_), "sample_ids")
    expect_error_message(open_reader(1), "sample_ids")
  })
  run_named_test("invalid queries and unsupported ploidy fail clearly", {
    expect_error_message(readPgenWindow(reader, "1", 200, 300), "No variants")
    expect_error_message(readPgenWindow(reader, "3", 100, 110), "No variants")
    expect_error_message(readPgenWindow(reader, "1", 110, 100), "bounds")
    for (start in list(NA_real_, Inf, 0, 1.5, "100", numeric())) {
      expect_error_message(readPgenWindow(reader, "1", start, 110), "bounds")
    }
    expect_error_message(readPgenWindow(reader, NA_character_, 100, 110), "chrom")
    expect_error_message(readPgenWindow(reader, character(), 100, 110), "chrom")
    for (chrom in c("X", "chrY", "MT", "M", "XY", as.character(23:26),
                    "contig1", "1.5", "0")) {
      expect_error_message(readPgenWindow(reader, chrom, 100, 110), "ploidy")
    }
  })
  run_named_test("local input checks reject cloud URIs and directories", {
    for (i in seq_along(paths)) {
      bad <- as.list(paths)
      bad[[i]] <- "gs://bucket/input"
      expect_error_message(do.call(openPgenReader, bad), "localization")
      bad[[i]] <- "s3://bucket/input"
      expect_error_message(do.call(openPgenReader, bad), "localization")
      bad[[i]] <- fixture_dir
      expect_error_message(do.call(openPgenReader, bad), "readable local")
      bad[[i]] <- file.path(fixture_dir, "absent")
      expect_error_message(do.call(openPgenReader, bad), "readable local")
    }
    Sys.chmod(paths[[3L]], "0000")
    tryCatch(expect_error_message(open_reader(), "readable local"),
             finally = Sys.chmod(paths[[3L]], "0600"))
  })
  run_named_test("Zstandard PVAR paths fail with a format-specific error", {
    zstd_path <- file.path(fixture_dir, "unsupported.pvar.zst")
    # Valid plain text at this path proves rejection occurs before file parsing.
    writeLines(pvar_lines, zstd_path)
    expect_error_message(openPgenReader(paths[[1L]], zstd_path, paths[[3L]]),
                         "Zstandard PVAR")
  })
  run_named_test("duplicate IIDs and metadata count mismatches fail", {
    writeLines(c(psam_lines, psam_lines[[2L]]), paths[[3L]])
    expect_error_message(open_reader(), "duplicate")
    writeLines(psam_lines[-length(psam_lines)], paths[[3L]])
    expect_error_message(open_reader(), "sample count")
    writeLines(c(psam_lines, "0\textra\t0"), paths[[3L]])
    expect_error_message(open_reader(), "sample count")
    writeLines(psam_lines, paths[[3L]])
    writeLines(pvar_lines[-length(pvar_lines)], paths[[2L]])
    expect_error_message(open_reader(), "variant count")
    writeLines(pvar_lines, paths[[2L]])
  })
  run_named_test("metadata is parsed once and reads reuse the open handle", {
    unlink(paths[2:3])
    for (i in seq_len(3L)) {
      expect_identical_value(readPgenWindow(reader, "1", 100, 110), expected)
      readPgenWindow(reader, "2", 90, 90)
    }
    writeLines(pvar_lines, paths[[2L]])
    writeLines(psam_lines, paths[[3L]])
  })
  run_named_test("invalid metadata fails before genotype reads", {
    invalid_pvar <- pvar_lines
    invalid_pvar[[7L]] <- "1\t0\tinvalid\tA\tG"
    writeLines(invalid_pvar, paths[[2L]])
    expect_error_message(open_reader(), "positions")
    invalid_pvar[[7L]] <- "chr\t102\tinvalid\tA\tG"
    writeLines(invalid_pvar, paths[[2L]])
    expect_error_message(open_reader(), "chromosomes")
    invalid_pvar <- pvar_lines
    invalid_pvar[[3L]] <- "#CHROM\tPOS\tID\tREF\tREF"
    writeLines(invalid_pvar, paths[[2L]])
    expect_error_message(open_reader(), "invalid columns")
    writeLines(pvar_lines, paths[[2L]])
    writeLines(c("#FID\tSAMPLE\tSEX", psam_lines[-1L]), paths[[3L]])
    expect_error_message(open_reader(), "invalid columns")
    writeLines(psam_lines, paths[[3L]])
  })
  run_named_test("sex chromosome records outside an autosomal query are allowed", {
    mixed_pvar <- pvar_lines
    mixed_pvar[[4L]] <- "Y\t90\tsex_chromosome\tA\tG"
    writeLines(mixed_pvar, paths[[2L]])
    mixed_reader <- open_reader()
    tryCatch({
      expect_identical_value(readPgenWindow(mixed_reader, "1", 100, 110), expected)
      expect_error_message(readPgenWindow(mixed_reader, "Y", 90, 90), "ploidy")
    }, finally = closePgenReader(mixed_reader))
    writeLines(pvar_lines, paths[[2L]])
  })
  run_named_test("multiallelic records fail only inside the query", {
    multi <- file.path(fixture_dir, c("multi.pgen", "multi.pvar", "multi.psam"))
    multi_writer <- pgenlibr::NewPgenWriter(multi[[1L]], 6L, 2L,
                                           nonref_flags = FALSE, allele_ct_limit = 3L)
    tryCatch({
      alleles <- matrix(c(0L, 0L, 0L, 1L, 1L, 2L, 2L, 2L,
                          NA_integer_, NA_integer_, 0L, 2L), nrow = 2L)
      pgenlibr::AppendAlleles(multi_writer, alleles, allele_ct = 3L)
      pgenlibr::AppendBiallelic(multi_writer, c(0L, 1L, 2L, 0L, 1L, 2L))
    }, finally = pgenlibr::ClosePgenWriter(multi_writer))
    writeLines(c("#CHROM\tPOS\tID\tREF\tALT", "1\t99\tmulti\tA\tC,G",
                 "1\t100\tbiallelic\tA\tG"), multi[[2L]])
    writeLines(psam_lines, multi[[3L]])
    multi_reader <- do.call(openPgenReader, as.list(multi))
    tryCatch({
      expect_error_message(readPgenWindow(multi_reader, "1", 99, 100), "multiallelic")
      want <- matrix(c(0L, 1L, 2L, 0L, 1L, 2L), nrow = 1L,
                     dimnames = list("chr1_100_A_G", sample_ids))
      expect_identical_value(readPgenWindow(multi_reader, "1", 100, 100), want)
    }, finally = closePgenReader(multi_reader))
  })
  run_named_test("close is idempotent and blocks subsequent reads", {
    pgen_handle <- reader$pgen
    pvar_handle <- reader$pvar
    closePgenReader(reader)
    closePgenReader(reader)
    expect_error_message(readPgenWindow(reader, "1", 100, 110), "closed")
    expect_error_message(pgenlibr::GetVariantCt(pgen_handle), "closed")
    expect_error_message(pgenlibr::GetVariantChrom(pvar_handle, 1L), "closed")
    expect_error_message(readPgenWindow(NULL, "1", 100, 110), "PgenReader")
    expect_error_message(closePgenReader(NULL), "PgenReader")
  })
  invisible(TRUE)
}
