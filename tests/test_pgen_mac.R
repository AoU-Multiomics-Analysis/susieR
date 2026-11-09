# Test MAC filtering with the real batch runner and PLINK. Model fitting is not
# needed here: stop at LoadData, after extraction and the PGEN reader checks.
benchmark_library <- Sys.getenv("BENCH_LIBRARY", unset = "")
if (nzchar(benchmark_library)) .libPaths(c(benchmark_library, .libPaths()))
args <- commandArgs(trailingOnly = TRUE)
repo_dir <- normalizePath(if (length(args)) args[[1L]] else getwd(), mustWork = TRUE)
scripts_dir <- Sys.getenv("SUSIER_SCRIPTS_PATH", file.path(repo_dir, "R", "scripts"))
runner <- new.env(parent = globalenv())
sys.source(file.path(scripts_dir, "run_pgen_susie_batch.R"), envir = runner)

run_mac_extraction_checks <- function() {
  previous_options <- options(warning.length = 8170)
  on.exit(options(previous_options), add = TRUE)
  plink <- Sys.getenv("PLINK2", unset = Sys.which("plink2"))
  stopifnot("Set PLINK2 to a real PLINK 2 executable" = nzchar(plink))
  directory <- tempfile("pgen MAC 'quoted' $literal ", tmpdir = dirname(tempdir()))
  dir.create(directory)
  passed <- FALSE
  on.exit(if (passed) unlink(directory, recursive = TRUE) else
    message("MAC fixture and PLINK logs retained at: ", directory), add = TRUE)

  samples <- as.character(101:108)
  variants <- c("monomorphic", "below_alt", "at_alt", "above_alt", "below_ref",
                "at_ref", "above_ref", "at_with_missing", "children_only", "outside_window")
  positions <- c(seq.int(950L, 1030L, by = 10L), 2000L)
  below <- c(2L, 1L, 1L, 0L, 0L, 0L, 0L, 0L)
  at <- c(2L, 1L, 1L, 1L, 0L, 0L, 0L, 0L)
  above <- c(2L, 2L, 1L, 1L, 0L, 0L, 0L, 0L)
  hardcalls <- rbind(rep(0L, 8L), below, at, above, 2L - below, 2L - at, 2L - above,
    c(2L, 2L, 1L, 0L, 0L, 0L, 0L, NA_integer_),
    c(0L, 0L, 0L, 0L, 2L, 1L, 1L, 1L), at)
  dimnames(hardcalls) <- list(variants, samples)
  # Hand-checked counts are 0,4,5,6,4,5,6,5,5,5. The final record is outside
  # the cis window. The four children have no phenotype columns.
  pgen <- file.path(directory, "input.pgen")
  pvar <- file.path(directory, "input.pvar")
  psam <- file.path(directory, "input.psam")
  writer <- pgenlibr::NewPgenWriter(pgen, length(samples), nrow(hardcalls), nonref_flags = FALSE)
  for (i in seq_len(nrow(hardcalls))) pgenlibr::AppendBiallelic(writer, hardcalls[i, ])
  pgenlibr::ClosePgenWriter(writer)
  readr::write_tsv(tibble::tibble(`#CHROM` = "1", POS = positions, ID = variants,
                                 REF = "A", ALT = "G"), pvar)
  readr::write_tsv(tibble::tibble(`#FID` = "F", IID = samples,
    PAT = c(rep("0", 4L), rep("101", 4L)),
    MAT = c(rep("0", 4L), rep("102", 4L)), SEX = 0L), psam)
  bed <- file.path(directory, "phenotypes.bed")
  writeLines(c("#chr\tstart\tend\tphenotype_id\t101\t102\t103\t104",
               "1\t999\t1001\ttrait1\t1\t2\t3\t4"), bed)
  qtl <- file.path(directory, "qtl.tsv")
  writeLines(c("phenotype_id\tqval", "trait1\t0.01"), qtl)
  covariates <- file.path(directory, "covariates.tsv")
  writeLines(c("ID\t101\t102\t103\t104", "pc1\t0\t1\t0\t1"), covariates)
  worklist <- file.path(directory, "batch.tsv")
  readr::write_tsv(tibble::tibble(batch_id = "batch_000001", chromosome = "1",
    gene_id = "gene1", outcome_id = "trait1", phenotype_id = "trait1", phenotype_pos = 1000L,
    requested_window_start = 900L, requested_window_end = 1100L,
    window_start = 900L, window_end = 1100L, region = "chr1:900-1100"), worklist)
  base_args <- c("--pgen", pgen, "--pvar", pvar, "--psam", psam, "--worklist", worklist,
    "--phenotype-bed", bed, "--qtl-results", qtl, "--covariates", covariates,
    "--missing-code=NA", "--cis-distance", "100", "--plink2", plink)
  parser <- optparse::OptionParser(option_list = runner$batch_options)
  runner$LoadData <- function(...) {
    stop(structure(list(message = "Extraction checked before phenotype loading", call = NULL),
      class = c("mac_extraction_complete", "error", "condition")))
  }
  extract <- function(label, mac = NULL, cli = FALSE, expected_error = NULL) {
    output <- file.path(directory, label)
    arguments <- c(base_args, "--output-dir", output,
      if (cli && !is.null(mac)) c("--mac", as.character(mac)))
    opt <- optparse::parse_args(parser, args = arguments, print_help_and_exit = FALSE)
    if (!cli && !is.null(mac)) opt$mac <- mac
    result <- tryCatch(runner$runPgenSusieBatch(opt), error = identity)
    if (!is.null(expected_error)) {
      stopifnot(inherits(result, "error"), !inherits(result, "mac_extraction_complete"))
      if (!grepl(expected_error, conditionMessage(result), fixed = TRUE)) {
        stop("Unexpected MAC extraction error: ", conditionMessage(result), call. = FALSE)
      }
      return(invisible(NULL))
    }
    if (!inherits(result, "mac_extraction_complete")) {
      stop("MAC batch failed before the extraction check: ", conditionMessage(result), call. = FALSE)
    }
    reader <- openPgenReader(paste0(file.path(output, "batch.genotypes"), ".pgen"),
      paste0(file.path(output, "batch.genotypes"), ".pvar"),
      paste0(file.path(output, "batch.genotypes"), ".psam"))
    on.exit(closePgenReader(reader), add = TRUE)
    stopifnot(identical(reader$sample_ids, samples),
      identical(reader$raw_sample_ids, samples), reader$raw_sample_count == 8L)
    list(ids = reader$variants$original_id,
      calls = readPgenWindow(reader, "1", 900L, 1100L),
      pvar = readLines(paste0(file.path(output, "batch.genotypes"), ".pvar")),
      psam = readLines(paste0(file.path(output, "batch.genotypes"), ".psam")))
  }

  # A missing --mac flag in the production PLINK call retains below-threshold
  # variants and fails here. These include both minor ALT and minor REF alleles.
  enabled <- extract("enabled", 5L)
  stopifnot("MAC 5 must keep counts equal to or above 5 across all PGEN samples" =
    identical(enabled$ids, c("at_alt", "above_alt", "at_ref", "above_ref", "at_with_missing", "children_only")))
  expected_rows <- c(3L, 4L, 6L, 7L, 8L, 9L)
  stopifnot(identical(unname(enabled$calls), unname(hardcalls[expected_rows, , drop = FALSE])))
  cli_enabled <- extract("cli_enabled", 5L, cli = TRUE)
  stopifnot(identical(cli_enabled, enabled))

  omitted <- extract("omitted")
  zero <- extract("zero", 0L, cli = TRUE)
  stopifnot("Omitted MAC and --mac 0 must preserve the same extracted data" = identical(omitted, zero),
    identical(omitted$ids, variants[1:9]),
    identical(unname(omitted$calls), unname(hardcalls[1:9, , drop = FALSE])))
  extract("all_removed", 17L, cli = TRUE, expected_error = "No variants remaining")
  passed <- TRUE
  message("Real PLINK MAC extraction checks passed")
}

run_mac_extraction_checks()
