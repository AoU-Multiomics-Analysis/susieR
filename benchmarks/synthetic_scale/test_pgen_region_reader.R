# Small, real PGEN checks for the benchmark reader. No production files are used.
# Source pgen_region_reader.R, then call run_pgen_reader_checks(tmpdir).
run_pgen_reader_checks <- function(tmpdir = tempdir()) {
  if (!exists("readPgenRegionHardcalls", mode = "function")) {
    stop("PGEN region reader is not defined", call. = FALSE)
  }
  if (!requireNamespace("pgenlibr", quietly = TRUE)) {
    stop("The benchmark checks require pgenlibr", call. = FALSE)
  }

  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
  fixture_dir <- tempfile("pgen-reader-checks-", tmpdir = tmpdir)
  dir.create(fixture_dir)
  on.exit(unlink(fixture_dir, recursive = TRUE), add = TRUE)
  paths <- file.path(fixture_dir, c("small.pgen", "small.pvar", "small.psam"))
  names(paths) <- c("pgen", "pvar", "psam")
  sample_ids <- c("s3", "001", "NA", "s2", "s1", "s4")
  # Rows before the query prevent region-relative indexes from passing.
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
  writer <- pgenlibr::NewPgenWriter(paths[["pgen"]], 6L, 12L,
                                   nonref_flags = FALSE)
  writer_open <- TRUE
  on.exit(if (writer_open) pgenlibr::ClosePgenWriter(writer), add = TRUE)
  for (i in seq_len(nrow(hardcalls))) {
    pgenlibr::AppendBiallelic(writer, hardcalls[i, ])
  }
  pgenlibr::ClosePgenWriter(writer)
  writer_open <- FALSE

  pvar_lines <- c(
    "##fileformat=PVARv1.2", "##fixture=benchmark-reader",
    "#CHROM\tPOS\tID\tREF\tALT",
    "2\t90\tother_chromosome\tA\tG", "1\t90\tbefore_query\tA\tG",
    "chr1\t100\tleft_boundary\tA\tG", "1\t101\tzero\tA\tG",
    "CHR1\t102\tzero_missing\tA\tG", "1\t103\tone_missing\tA\tG",
    "chr1\t104\ttwo_missing\tA\tG", "1\t105\tall_missing\tA\tG",
    "chr1\t106\tmixed\tA\tG", "1\t107\tNA\tNA\tG",
    "1\t110\tright_boundary\tA\tG", "1\t111\tafter_query\tA\tG"
  )
  psam_lines <- c("#FID\tIID\tSEX", paste("0", sample_ids, "0", sep = "\t"))
  writeLines(pvar_lines, paths[["pvar"]])
  writeLines(psam_lines, paths[["psam"]])
  read_region <- function(...) {
    readPgenRegionHardcalls(paths[["pgen"]], paths[["pvar"]], paths[["psam"]], ...)
  }
  checks <- list()
  check <- function(name, passed) {
    if (!isTRUE(passed)) stop("PGEN reader check failed: ", name, call. = FALSE)
    checks[[length(checks) + 1L]] <<- tibble::tibble(check = name, passed = TRUE)
  }
  expect_error <- function(name, expression, pattern) {
    result <- tryCatch(force(expression), error = identity)
    check(name, inherits(result, "error") && grepl(pattern, conditionMessage(result),
                                                   ignore.case = TRUE))
  }

  expected_ids <- c("chr1_100_A_G", "chr1_101_A_G", "chr1_102_A_G",
                    "chr1_103_A_G", "chr1_104_A_G", "chr1_105_A_G",
                    "chr1_106_A_G", "chr1_107_NA_G", "chr1_110_A_G")
  expected <- hardcalls[3:11, , drop = FALSE]
  dimnames(expected) <- list(expected_ids, sample_ids)
  actual <- read_region("1", 100, 110, as_double = FALSE)
  check("global indexes, inclusive bounds, orientation and exact hardcalls",
        identical(actual, expected))
  check("missing hardcalls remain missing", identical(is.na(actual), is.na(expected)))
  check("metadata strings preserve IID NA, leading zeros and allele NA",
        identical(colnames(actual)[2:3], c("001", "NA")) &&
          "chr1_107_NA_G" %in% rownames(actual))
  double_expected <- expected
  storage.mode(double_expected) <- "double"
  check("default result is a plain double matrix",
        identical(read_region("chr1", 100, 110), double_expected))
  check("chromosome prefixes are case insensitive",
        identical(read_region("CHRchr1", 100, 110), double_expected))
  reordered_ids <- c("s1", "001", "NA")
  check("sample subset is restored to requested order",
        identical(read_region("1", 100, 110, reordered_ids, FALSE),
                  expected[, reordered_ids, drop = FALSE]))
  check("single variant keeps both matrix dimensions",
        identical(read_region("1", 100, 100, "s1", FALSE),
                  expected[1L, "s1", drop = FALSE]))
  expect_error("empty region fails", read_region("1", 200, 300), "no variants")
  expect_error("unknown sample fails", read_region("1", 100, 110, "absent"), "sample")
  expect_error("duplicate requested samples fail",
               read_region("1", 100, 110, c("s1", "s1")), "duplicate")
  expect_error("empty sample subset fails", read_region("1", 100, 110, character()),
               "sample")
  expect_error("reversed bounds fail", read_region("1", 110, 100), "start|end|bound")

  writeLines(c(psam_lines, psam_lines[[2L]]), paths[["psam"]])
  expect_error("duplicate PSAM IID fails", read_region("1", 100, 110), "duplicate")
  writeLines(psam_lines[-length(psam_lines)], paths[["psam"]])
  expect_error("PSAM and PGEN sample counts must match", read_region("1", 100, 110),
               "sample.*count|count.*sample")
  writeLines(psam_lines, paths[["psam"]])

  # A true multiallelic record must fail if it is inside the query only.
  multi_paths <- file.path(fixture_dir, c("multi.pgen", "multi.pvar", "multi.psam"))
  multi_writer <- pgenlibr::NewPgenWriter(multi_paths[[1L]], 6L, 2L,
                                         nonref_flags = FALSE, allele_ct_limit = 3L)
  multi_open <- TRUE
  on.exit(if (multi_open) pgenlibr::ClosePgenWriter(multi_writer), add = TRUE)
  allele_codes <- matrix(c(0L, 0L, 0L, 1L, 1L, 2L, 2L, 2L,
                           NA_integer_, NA_integer_, 0L, 2L), nrow = 2L)
  pgenlibr::AppendAlleles(multi_writer, allele_codes, allele_ct = 3L)
  pgenlibr::AppendBiallelic(multi_writer, c(0L, 1L, 2L, 0L, 1L, 2L))
  pgenlibr::ClosePgenWriter(multi_writer)
  multi_open <- FALSE
  writeLines(c("#CHROM\tPOS\tID\tREF\tALT", "1\t99\tmulti\tA\tC,G",
               "1\t100\tbiallelic\tA\tG"), multi_paths[[2L]])
  writeLines(psam_lines, multi_paths[[3L]])
  read_multi <- function(start, end) {
    readPgenRegionHardcalls(multi_paths[[1L]], multi_paths[[2L]], multi_paths[[3L]],
                           "1", start, end, as_double = FALSE)
  }
  expect_error("query multiallelic record fails", read_multi(99, 100), "multiallelic")
  multi_expected <- matrix(c(0L, 1L, 2L, 0L, 1L, 2L), nrow = 1L,
                           dimnames = list("chr1_100_A_G", sample_ids))
  check("multiallelic record outside query does not remove biallelic records",
        identical(read_multi(100, 100), multi_expected))

  # The old filter uses -1 in its sum, before imputation. This is observable.
  legacy <- expected
  legacy[3L, 2L] <- -1L
  legacy[4L, 2L] <- -1L
  legacy[6L, c(1L, 3L, 5L)] <- -1L
  nonzero <- function(x) rowSums(round(x, 0), na.rm = TRUE) != 0
  check("NA filter boundary expectations",
        identical(unname(nonzero(expected)),
                  c(TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE)))
  check("mixed legacy missing codes change three filter decisions",
        identical(unname(nonzero(legacy)),
                  c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE)) &&
          sum(nonzero(legacy) != nonzero(expected)) == 3L)
  legacy[legacy == -1L] <- NA_integer_
  check("normalization before filtering restores matched TSV and PGEN values",
        identical(legacy, actual) && identical(nonzero(legacy), nonzero(actual)))

  dplyr::bind_rows(checks)
}
