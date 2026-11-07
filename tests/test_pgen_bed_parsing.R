# Metadata selection and row splitting must use the same BED field rules.
run_pgen_bed_parsing_tests <- function(repo_dir) {
  suppressPackageStartupMessages(library(dplyr))
  functions_dir <- Sys.getenv("SUSIER_FUNCTIONS_PATH",
                              unset = file.path(repo_dir, "R", "utils"))
  source(file.path(functions_dir, "ImportFunctions.R"))
  source(file.path(functions_dir, "PgenBatchFunctions.R"))
  root <- tempfile("pgen BED parsing ")
  dir.create(root)
  old <- setwd(root)
  on.exit({ setwd(old); unlink(root, recursive = TRUE) }, add = TRUE)
  error <- function(expr, pattern, label) {
    result <- tryCatch({ force(expr); NULL }, error = identity)
    if (!inherits(result, "error") ||
        !grepl(pattern, conditionMessage(result), ignore.case = TRUE)) {
      stop(label, ": expected ", pattern, "; got ",
           if (inherits(result, "error")) conditionMessage(result) else "no error",
           call. = FALSE)
    }
    invisible(result)
  }
  write_bytes <- function(lines, path, eol = "\n") {
    writeBin(charToRaw(paste0(paste(lines, collapse = eol), eol)), path)
  }
  header <- "#chr\tstart\tend\tphenotype_id\tS 2\tS1"
  rows <- c(
    'chr1\t10\t12\t"p_quote"\t1.00000000\tNA',
    "chr1\t20\t22\t p_left\t1e-10\t2",
    "chr1\t30\t32\tp_right \t3\t4",
    'chr1\t40\t42\t"p""quote"\t5\t6',
    "chr1\t50\t52\tliteral#id\t7\t8",
    "chr1\t60\t62\t0007\t9\t10",
    "chr1\t70\t72\tp_drop\t11\t12")
  # read.table() accepts single-quoted QTL IDs containing a literal double quote
  # or a comment marker. This fixture isolates BED parsing from QTL parsing.
  writeLines(c("phenotype_id\tqval", "p_quote\t0.01", "p_left\t0.01",
               "p_right\t0.01", "'p\"quote'\t0.01", "'literal#id'\t0.01",
               "0007\t0.01", "p_drop\t0.9"), "qtl.tsv")
  expected_ids <- c("p_quote", "p_left", "p_right", 'p"quote', "literal#id", "0007")
  check_plan <- function(input, output, expected_header = header) {
    planPgenBatches(input, "qtl.tsv", "chr1", 5, output, batch_size = 2)
    expected <- readr::read_tsv(file.path(output, "expected_outcomes.tsv"),
      col_types = readr::cols(.default = readr::col_character()),
      show_col_types = FALSE)
    stopifnot(identical(expected$phenotype_id, expected_ids),
              identical(expected$batch_id, c("batch_000001", "batch_000001",
                "batch_000002", "batch_000002", "batch_000003", "batch_000003")),
              identical(expected$phenotype_pos, c("11", "21", "31", "41", "51", "61")))
    paths <- readLines(file.path(output, "batch_phenotype_beds.txt"), warn = FALSE)
    stopifnot(length(paths) == 3L,
              identical(readLines(paths[[1L]], warn = FALSE), c(expected_header, rows[1:2])),
              identical(readLines(paths[[2L]], warn = FALSE), c(expected_header, rows[3:4])),
              identical(readLines(paths[[3L]], warn = FALSE), c(expected_header, rows[5:6])))
  }
  # Reintroducing raw-ID matching must fail this quoted/padded-ID case.
  write_bytes(c(header, rows), "quoted.bed")
  check_plan("quoted.bed", "quoted plan")

  # Compression must not select different field rules or depend on the suffix.
  connection <- gzfile("quoted.bed.gz", "wb")
  writeBin(readBin("quoted.bed", "raw", n = file.info("quoted.bed")$size), connection)
  close(connection)
  check_plan("quoted.bed.gz", "gzip plan")
  for (suffix in c("bgz", "bgzf", "blob")) {
    input <- paste0("quoted.bed.", suffix)
    Rsamtools::bgzip("quoted.bed", dest = input, overwrite = TRUE)
    check_plan(input, paste0(suffix, " plan"))
  }

  # Line endings and preambles must not change the header or selected row text.
  write_bytes(c(header, rows), "crlf.bed", eol = "\r\n")
  check_plan("crlf.bed", "crlf plan")
  bom_header <- paste0("\ufeff", header)
  write_bytes(c(bom_header, rows), "bom.bed")
  # readLines() removes the initial UTF-8 BOM before returning header text.
  check_plan("bom.bed", "bom plan")
  preamble <- c("", "# BED comment", 'track name="phenotypes"',
                "browser position chr1:1-100", "", header, rows)
  write_bytes(preamble, "preamble.bed")
  check_plan("preamble.bed", "preamble plan")
  Rsamtools::bgzip("preamble.bed", dest = "preamble.bed.bgz", overwrite = TRUE)
  check_plan("preamble.bed.bgz", "compressed preamble plan")
  write_bytes(c(header, rows[1:3], "", rows[4:7]), "blank data row.bed")
  check_plan("blank data row.bed", "blank data row plan")

  # Selected rows on both sides of the streaming boundary must appear once.
  many_ids <- paste0("m", seq_len(1201L))
  many_rows <- paste("chr1", seq_len(1201L), seq_len(1201L) + 1L,
                     many_ids, "1.00000000", "NA", sep = "\t")
  many_rows[[1000L]] <- 'chr1\t1000\t1001\t"m1000"\t1.00000000\tNA'
  many_rows[[1001L]] <- "chr1\t1001\t1002\t m1001 \t1.00000000\tNA"
  write_bytes(c(header, many_rows), "boundary.bed")
  readr::write_tsv(tibble::tibble(phenotype_id = c("m1000", "m1001", "m1201"),
                                 qval = .01), "boundary.tsv")
  planPgenBatches("boundary.bed", "boundary.tsv", "chr1", 5, "boundary plan")
  path <- readLines("boundary plan/batch_phenotype_beds.txt", warn = FALSE)
  stopifnot(length(path) == 1L,
            identical(readLines(path, warn = FALSE),
                      c(header, many_rows[c(1000L, 1001L, 1201L)])))

  # Invalid delimiters and row widths must fail before output BEDs are made.
  malformed <- list(
    spaces = gsub("\t", " ", c(header, rows[[1L]]), fixed = TRUE),
    short_header = c("#chr\tstart\tend", rows[[1L]]),
    wrong_header = c("#chr\treference\tbuild\tphenotype_id\tS 2\tS1", rows[[1L]]),
    short_row = c(header, 'chr1\t10\t12\t"p_quote"\t1'),
    long_row = c(header, paste0(rows[[1L]], "\textra")),
    missing_id = c(header, "chr1\t10\t12\t\t1\t2"),
    unbalanced_quote = c(header, 'chr1\t10\t12\t"p_quote\t1\t2'),
    text_after_quote = c(header, 'chr1\t10\t12\t"p_quote"extra\t1\t2'))
  for (name in names(malformed)) {
    input <- paste0(name, ".bed")
    output <- paste0(name, " plan")
    write_bytes(malformed[[name]], input)
    pattern <- if (name %in% c("unbalanced_quote", "text_after_quote")) {
      "quote|metadata|row"
    } else "tab|column|header|row|nonempty"
    error(planPgenBatches(input, "qtl.tsv", "chr1", 5, output),
          pattern, name)
    stopifnot(!file.exists(file.path(output, "batch_phenotype_beds.txt")))
  }
  write_bytes(c(header, rows[[1L]], 'chr1\t20\t22\t p_quote \t3\t4'), "duplicate.bed")
  error(planPgenBatches("duplicate.bed", "qtl.tsv", "chr1", 5, "duplicate plan"),
        "duplicate", "IDs that become equal after decoding")

  # Change the real input between passes. Missing selected rows must still fail.
  write_bytes(c(header, "chr1\t10\t12\tp1\t1\t2", "chr1\t20\t22\tp2\t3\t4"),
              "changed.bed")
  readr::write_tsv(tibble::tibble(phenotype_id = c("p1", "p2"), qval = .01), "changed.tsv")
  changed_plan <- planPgenBatches
  changed_environment <- new.env(parent = environment(changed_plan))
  import_qtl <- importQtlmapPermutedPvalues
  changed_environment$importQtlmapPermutedPvalues <- function(path) {
    result <- import_qtl(path)
    write_bytes(c(header, "chr1\t20\t22\tp2\t3\t4"), "changed.bed")
    result
  }
  environment(changed_plan) <- changed_environment
  error(changed_plan("changed.bed", "changed.tsv", "chr1", 5, "changed plan"),
        "changed|missing|split", "Selected row removed between passes")
  # The same IDs and row count do not prove that planned positions are current.
  write_bytes(c(header, "chr1\t10\t12\tp1\t1\t2", "chr1\t20\t22\tp2\t3\t4"),
              "changed.bed")
  changed_environment$importQtlmapPermutedPvalues <- function(path) {
    result <- import_qtl(path)
    write_bytes(c(header, "chr1\t100\t102\tp1\t1\t2", "chr1\t20\t22\tp2\t3\t4"),
                "changed.bed")
    result
  }
  error(changed_plan("changed.bed", "changed.tsv", "chr1", 5, "changed position plan"),
        "changed|metadata|position", "Selected position changed between passes")
  message("PGEN BED parsing tests passed: matching IDs, exact rows, compression, and input checks")
  invisible(TRUE)
}
