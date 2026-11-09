# Real PLINK fixtures test shared QC before any phenotype sample selection.
run_pgen_qc_tests <- function(repo_dir) {
  source(file.path(repo_dir, "tests", "test_helpers.R"), local = TRUE)
  functions_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH",
                               unset = file.path(repo_dir, "R", "utils"))
  qc_path <- file.path(functions_path, "PgenQcFunctions.R")
  if (!file.exists(qc_path)) stop("Missing shared PGEN QC implementation", call. = FALSE)
  source(file.path(functions_path, "PgenFunctions.R"), local = TRUE)
  source(qc_path, local = TRUE)
  plink <- Sys.getenv("PLINK2", unset = Sys.which("plink2"))
  if (!nzchar(plink)) stop("PLINK2 must name a real PLINK 2 executable", call. = FALSE)
  fixture_dir <- tempfile("shared pgen qc '")
  dir.create(fixture_dir)
  on.exit(unlink(fixture_dir, recursive = TRUE), add = TRUE)
  run_plink <- function(args) {
    output <- system2(plink, vapply(args, shQuote, character(1)), stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(output, "status"))) stop(paste(output, collapse = "\n"), call. = FALSE)
    invisible(output)
  }
  sample_ids <- paste0("s", seq_len(20L))
  zero <- rep("0|0:0", 20L)
  genotypes <- rep(list(zero), 11L)
  genotypes[[2L]][c(1L, 20L)] <- c("0|1:1", "./.:.")
  genotypes[[3L]][19:20] <- "./.:."
  genotypes[[4L]][1:6] <- "0|1:1"
  genotypes[[5L]][[1L]] <- "0|1:1"
  genotypes[[6L]] <- rep("1|1:2", 20L)
  genotypes[[7L]][[1L]] <- "./.:0.4"
  genotypes[[8L]][[20L]] <- "1|0:1"
  genotypes[[9L]] <- c("0|1:1", rep("1|1:2", 19L))
  genotypes[[10L]] <- rep("./.:.", 20L)
  # Multiallelic records carry real allele indices; do not import DS for them.
  genotypes[[11L]] <- rep("0/2:.", 20L)
  raw_ids <- c("dup", "dup", ".", ".", "long", "stale", "dosage", "child", "minor_ref", "sex", "multi")
  alleles <- c(rep("C", 4L), paste(rep("A", 400L), collapse = ""), rep("C", 5L), "C,T")
  vcf <- file.path(fixture_dir, "input.vcf")
  header <- c("##fileformat=VCFv4.2", "##contig=<ID=1>", "##contig=<ID=X>",
              "##INFO=<ID=MAC,Number=1,Type=Integer,Description=\"Stale metadata\">",
              "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
              "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage\">",
              paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_ids), collapse = "\t"))
  records <- vapply(seq_along(genotypes), function(index) {
    paste(c("1", index, raw_ids[[index]], "G", alleles[[index]],
            ".", "PASS", if (index %in% c(1L, 6L)) "MAC=999" else ".", "GT:DS",
            genotypes[[index]]), collapse = "\t")
  }, character(1))
  # Import the biallelic phased and dosage records, then append a multiallelic
  # hardcall record with pgenlibr in a separate preservation fixture below.
  writeLines(c(header, records[1:10]), vcf)
  prefix <- file.path(fixture_dir, "input")
  run_plink(c("--vcf", vcf, "dosage=DS", "--make-pgen", "--out", prefix,
              "--threads", "1", "--memory", "640"))
  pvar_lines <- readLines(paste0(prefix, ".pvar"))
  pvar_rows <- which(!startsWith(pvar_lines, "#"))
  pvar_lines[pvar_rows[[10L]]] <- sub("^1\t", "X\t", pvar_lines[pvar_rows[[10L]]])
  writeLines(pvar_lines, paste0(prefix, ".pvar"))
  psam <- c("#IID\tPAT\tMAT\tSEX", vapply(seq_along(sample_ids), function(index) {
    paste(sample_ids[[index]], if (index == 20L) "s1" else "0",
          if (index == 20L) "s2" else "0", "NA", sep = "\t")
  }, character(1)))
  writeLines(psam, paste0(prefix, ".psam"))
  opts <- function(name, mac = 0, min_call_rate = .95) {
    list(pgen = paste0(prefix, ".pgen"), pvar = paste0(prefix, ".pvar"),
         psam = paste0(prefix, ".psam"), mac = mac, min_call_rate = min_call_rate,
         output_dir = file.path(fixture_dir, name), plink2 = plink, plink_memory_mb = 640)
  }
  read_counts <- function(output) readr::read_tsv(output$counts, show_col_types = FALSE,
                                                na = "NA", progress = FALSE)
  run_named_test("QC thresholds reject invalid numeric values", {
    validatePgenQcThresholds(0, .95)
    for (value in list(-1, .5, NA_real_, NaN, Inf, "1", c(0, 1), NULL)) {
      expect_error_message(validatePgenQcThresholds(value, .95), "mac")
    }
    for (value in list(-1, 1.01, NA_real_, NaN, Inf, ".95", c(0, 1), NULL)) {
      expect_error_message(validatePgenQcThresholds(0, value), "min_call_rate")
    }
  })
  result <- filterPgenFile(opts("default qc"))
  run_named_test("hardcall AN and MAC ignore dosage and stale PVAR annotations", {
    counts <- read_counts(result)
    stopifnot(nrow(counts) == 10L,
              identical(counts$AN[1:9], c(40, 38, 36, 40, 40, 40, 38, 40, 40)),
              identical(counts$MAC[1:9], c(0, 1, 0, 6, 1, 0, 0, 1, 1)),
              counts$AC[[9L]] == 39, counts$MAC[[8L]] == 1,
              counts$CALL_RATE[[2L]] == .95, counts$KEEP[[2L]],
              !counts$KEEP[[3L]], counts$KEEP[[7L]],
              !counts$qc_applied[[10L]], counts$KEEP[[10L]], is.na(counts$AN[[10L]]),
              identical(counts$variant_id[1:4], c("dup", "dup", ".", ".")))
    report <- jsonlite::read_json(result$report, simplifyVector = TRUE)
    stopifnot(report$schema_version == 1, report$count_source == "pgen_hardcalls",
              report$sample_scope == "all_input_pgen_samples", report$include_nonfounders,
              report$sample_count == 20, report$variant_count_before == 10,
              report$variant_count_after == 9, report$thresholds$mac == 0,
              report$thresholds$min_call_rate == .95,
              identical(report$filtered_files$pgen$sha256,
                        digest::digest(file = result$pgen, algo = "sha256", serialize = FALSE)),
              identical(report$counts_file$sha256,
                        digest::digest(file = result$counts, algo = "sha256", serialize = FALSE)))
  })
  run_named_test("QC preserves raw PVAR lines, sample order, dosage and phase", {
    original_lines <- readLines(paste0(prefix, ".pvar"))
    metadata <- original_lines[startsWith(original_lines, "#")]
    data <- original_lines[!startsWith(original_lines, "#")]
    expect_identical_value(readLines(result$pvar), c(metadata, data[-3L]))
    expect_identical_value(readBin(result$psam, "raw", n = file.info(result$psam)$size),
                           readBin(paste0(prefix, ".psam"), "raw", n = file.info(paste0(prefix, ".psam"))$size))
    for (name in c("input", "filtered")) {
      pgen_prefix <- if (name == "input") prefix else sub("\\.pgen$", "", result$pgen)
      run_plink(c("--pfile", pgen_prefix, "--export", "vcf", "vcf-dosage=DS", "--out",
                  file.path(fixture_dir, paste0("export_", name)), "--threads", "1", "--memory", "640"))
    }
    original <- readLines(file.path(fixture_dir, "export_input.vcf"))
    filtered <- readLines(file.path(fixture_dir, "export_filtered.vcf"))
    expect_identical_value(filtered[!startsWith(filtered, "#")],
                           original[!startsWith(original, "#")][-3L])
  })
  run_named_test("QC accepts leading whitespace and preserves retained PVAR bytes", {
    leading_opt <- opts("leading whitespace")
    leading_opt$pvar <- file.path(fixture_dir, "leading.pvar")
    leading_lines <- readLines(paste0(prefix, ".pvar"))
    data_rows <- which(!startsWith(leading_lines, "#"))
    leading_lines[data_rows] <- paste0(rep(c("  ", "\t "), length.out = length(data_rows)),
                                      leading_lines[data_rows])
    writeLines(leading_lines, leading_opt$pvar)
    leading <- filterPgenFile(leading_opt)
    expected_pvar <- file.path(fixture_dir, "expected-leading.pvar")
    writeLines(leading_lines[-data_rows[[3L]]], expected_pvar)
    expect_identical_value(readBin(leading$pvar, "raw", n = file.info(leading$pvar)$size),
                           readBin(expected_pvar, "raw", n = file.info(expected_pvar)$size))
    stopifnot(identical(read_counts(leading)$KEEP, read_counts(result)$KEEP))
  })
  run_named_test("minimum MAC filters minor reference alleles and exact row IDs", {
    mac1 <- filterPgenFile(opts("mac1", 1))
    counts <- read_counts(mac1)
    expect_identical_value(which(counts$KEEP), c(2L, 4L, 5L, 8L, 9L, 10L))
    expect_identical_value(readLines(mac1$pvar)[!startsWith(readLines(mac1$pvar), "#")],
                           readLines(paste0(prefix, ".pvar"))[!startsWith(readLines(paste0(prefix, ".pvar")), "#")][c(2L, 4L, 5L, 8L, 9L, 10L)])
    expect_error_message(filterPgenFile(opts("all fail", 100)), "No supported variants remain")
  })
  run_named_test("QC retains multiallelic records outside supported scope", {
    multi_prefix <- file.path(fixture_dir, "multi")
    writer <- pgenlibr::NewPgenWriter(paste0(multi_prefix, ".pgen"), 20L, 2L,
                                     allele_ct_limit = 3L, nonref_flags = FALSE)
    tryCatch({
      pgenlibr::AppendAlleles(writer, matrix(rep(c(0L, 2L), 20L), nrow = 2L), allele_ct = 3L)
      pgenlibr::AppendBiallelic(writer, c(1L, rep(0L, 19L)))
    }, finally = pgenlibr::ClosePgenWriter(writer))
    writeLines(c("#CHROM\tPOS\tID\tREF\tALT", "1\t1\t.\tG\tC,T", "1\t2\t.\tG\tC"), paste0(multi_prefix, ".pvar"))
    file.copy(paste0(prefix, ".psam"), paste0(multi_prefix, ".psam"))
    multi_opt <- opts("multi qc", 1)
    for (extension in c("pgen", "pvar", "psam")) multi_opt[[extension]] <- paste0(multi_prefix, ".", extension)
    multi <- filterPgenFile(multi_opt)
    counts <- read_counts(multi)
    stopifnot(identical(counts$qc_applied, c(FALSE, TRUE)), all(counts$KEEP), is.na(counts$MAC[[1L]]))
    expect_identical_value(readLines(multi$pvar), readLines(paste0(multi_prefix, ".pvar")))
  })
  run_named_test("QC streams across chunk boundaries and counts all-missing rows", {
    stream_prefix <- file.path(fixture_dir, "stream")
    writer <- pgenlibr::NewPgenWriter(paste0(stream_prefix, ".pgen"), 20L, 1005L,
                                     nonref_flags = FALSE)
    tryCatch({
      for (index in seq_len(1005L)) {
        values <- if (index %% 2L == 0L) rep(NA_integer_, 20L) else c(1L, rep(0L, 19L))
        pgenlibr::AppendBiallelic(writer, values)
      }
    }, finally = pgenlibr::ClosePgenWriter(writer))
    writeLines(c("#CHROM\tPOS\tID\tREF\tALT",
                 paste("1", seq_len(1005L), ".", "G", "C", sep = "\t")), paste0(stream_prefix, ".pvar"))
    file.copy(paste0(prefix, ".psam"), paste0(stream_prefix, ".psam"))
    stream_opt <- opts("stream qc", 1)
    for (extension in c("pgen", "pvar", "psam")) stream_opt[[extension]] <- paste0(stream_prefix, ".", extension)
    stream <- filterPgenFile(stream_opt)
    counts <- read_counts(stream)
    stopifnot(nrow(counts) == 1005L, identical(which(counts$KEEP), seq.int(1L, 1005L, 2L)),
              all(counts$AN[seq.int(2L, 1004L, 2L)] == 0),
              all(counts$CALL_RATE[seq.int(2L, 1004L, 2L)] == 0))
    stream_opt$output_dir <- file.path(fixture_dir, "disabled thresholds")
    stream_opt$mac <- 0
    stream_opt$min_call_rate <- 0
    disabled <- filterPgenFile(stream_opt)
    stopifnot(all(read_counts(disabled)$KEEP))
  })
  run_named_test("QC validates localization and memory before PLINK", {
    invalid <- opts("uri")
    invalid$pgen <- "gs://bucket/input.pgen"
    expect_error_message(filterPgenFile(invalid), "localization error")
    invalid <- opts("missing")
    invalid$pvar <- file.path(fixture_dir, "absent.pvar")
    expect_error_message(filterPgenFile(invalid), "readable local file")
    invalid <- opts("memory")
    invalid$plink_memory_mb <- .5
    expect_error_message(filterPgenFile(invalid), "plink_memory_mb")
    invalid <- opts("threshold")
    invalid$mac <- .5
    expect_error_message(filterPgenFile(invalid), "mac")
  })
  run_named_test("filter CLI applies defaults and rejects invalid thresholds", {
    cli <- file.path(Sys.getenv("SUSIER_SCRIPTS_PATH",
                    unset = file.path(repo_dir, "R", "scripts")), "filter_pgen_susie.R")
    args <- c(cli, "--pgen", paste0(prefix, ".pgen"), "--pvar", paste0(prefix, ".pvar"),
              "--psam", paste0(prefix, ".psam"), "--output-dir", file.path(fixture_dir, "cli"),
              "--plink2", plink, "--plink-memory-mb", "640")
    output <- system2(file.path(R.home("bin"), "Rscript"), vapply(args, shQuote, character(1)), stdout = TRUE, stderr = TRUE,
                      env = paste0("SUSIER_FUNCTIONS_PATH=", shQuote(functions_path)))
    stopifnot(is.null(attr(output, "status")))
    report <- jsonlite::read_json(file.path(fixture_dir, "cli", "qc.json"), simplifyVector = TRUE)
    stopifnot(report$thresholds$mac == 0, report$thresholds$min_call_rate == .95)
    output <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
      vapply(c(args, "--mac", "0.5"), shQuote, character(1)), stdout = TRUE, stderr = TRUE,
      env = paste0("SUSIER_FUNCTIONS_PATH=", shQuote(functions_path))))
    stopifnot(!is.null(attr(output, "status")), any(grepl("mac", output, fixed = TRUE)))
  })
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  benchmark_library <- Sys.getenv("BENCH_LIBRARY", unset = "")
  if (nzchar(benchmark_library)) .libPaths(c(benchmark_library, .libPaths()))
  args <- commandArgs(trailingOnly = TRUE)
  run_pgen_qc_tests(normalizePath(if (length(args)) args[[1L]] else getwd()))
}
