# Real native PGEN preparation, shared cohort QC and report validation.
benchmark_library <- Sys.getenv("BENCH_LIBRARY", unset = "")
if (nzchar(benchmark_library)) .libPaths(c(benchmark_library, .libPaths()))
repo_dir <- normalizePath(getwd())
functions_dir <- Sys.getenv("SUSIER_FUNCTIONS_PATH", file.path(repo_dir, "R", "utils"))
scripts_dir <- Sys.getenv("SUSIER_SCRIPTS_PATH", file.path(repo_dir, "R", "scripts"))
source(file.path(repo_dir, "tests", "test_helpers.R"))
source(file.path(scripts_dir, "run_pgen_susie_batch.R"))
source(file.path(functions_dir, "PgenPreparationFunctions.R"))
run_preparation_tests <- function() {
  directory <- tempfile("pgen preparation 'quoted' $literal ")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  plink <- Sys.getenv("PLINK2", Sys.which("plink2"))
  stopifnot(nzchar(plink))
  prefix <- file.path(directory, "input")
  ids <- c("s6", "s3", "s1", "s5", "s2", "s4")
  positions <- c(90L, 100L, 110L, 150L, 195L, 200L, 205L, 300L, 100L)
  gt <- matrix(rep(c(0L, 1L, 2L, 0L, 1L, 2L), length(positions)), nrow = length(positions), byrow = TRUE)
  # The only minor allele of the final variant is in a nonfounder.
  gt[8L, ] <- c(0L, 0L, 0L, 0L, 0L, 1L)
  writer <- pgenlibr::NewPgenWriter(paste0(prefix, ".pgen"), length(ids), nrow(gt), nonref_flags = FALSE)
  for (i in seq_len(nrow(gt))) pgenlibr::AppendBiallelic(writer, gt[i, ])
  pgenlibr::ClosePgenWriter(writer)
  pvar_header <- c('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
  pvar_rows <- vapply(seq_along(positions), function(index) paste(c(
    if (index == 9L) "chr2" else "chr1", positions[[index]],
    if (index <= 2L) "duplicate" else paste0("v", index), "A", "G",
    "0.123456789", "PASS", "AF=0.123456789"), collapse = "\t"), character(1))
  writeLines(c(pvar_header, pvar_rows), paste0(prefix, ".pvar"))
  psam_lines <- c("#FID\tIID\tPAT\tMAT\tSEX\tprecise_phenotype", vapply(seq_along(ids), function(index)
    paste(c("0", ids[[index]], if (index == 6L) "s6" else "0", if (index == 6L) "s3" else "0",
      "NA", "0.000123456789"), collapse = "\t"), character(1)))
  writeLines(psam_lines, paste0(prefix, ".psam"))
  work <- function(id, chromosome, centers) tibble::tibble(batch_id = id, chromosome = chromosome,
    gene_id = paste0("g", centers), outcome_id = paste0("o", centers), phenotype_id = paste0("o", centers),
    phenotype_pos = centers, requested_window_start = centers - 10L, requested_window_end = centers + 10L,
    window_start = centers - 10L, window_end = centers + 10L,
    region = paste0(chromosome, ":", centers - 10L, "-", centers + 10L))
  worklists <- file.path(directory, c("first.tsv", "second.tsv", "third.tsv", "fourth.tsv"))
  readr::write_tsv(work("batch_000003", "chr1", c(100L, 105L, 200L)), worklists[[1L]])
  readr::write_tsv(work("batch_000008", "1", 300L), worklists[[2L]])
  readr::write_tsv(work("batch_000009", "1", 200L), worklists[[3L]])
  readr::write_tsv(work("batch_000010", "chr2", 100L), worklists[[4L]])
  lists <- file.path(directory, "worklists.txt")
  writeLines(worklists, lists)
  opt <- list(pgen = paste0(prefix, ".pgen"), pvar = paste0(prefix, ".pvar"), psam = paste0(prefix, ".psam"),
    worklists_file = lists, cis_distance = 10L, mac = 1L, min_call_rate = .95,
    plink2 = plink, plink_memory_mb = 640L, output_dir = file.path(directory, "prepared"))
  original_filter <- filterPgenFile
  calls <- 0L
  filterPgenFile <<- function(opt) { calls <<- calls + 1L; original_filter(opt) }
  on.exit(assign("filterPgenFile", original_filter, envir = .GlobalEnv), add = TRUE)
  stale <- opt
  stale$output_dir <- file.path(directory, "stale")
  dir.create(stale$output_dir)
  writeLines("Previous preparation output", file.path(stale$output_dir, "previous.txt"))
  stale_error <- tryCatch(preparePgenBatches(stale), error = identity)
  stopifnot("Nonempty preparation output directories must be rejected" = inherits(stale_error, "error"))
  expect_error_message(stop(stale_error), "output_dir must be empty")
  stopifnot(calls == 0L,
    identical(list.files(stale$output_dir), "previous.txt"))
  result <- preparePgenBatches(opt)
  run_named_test("one shared QC and one native shard per planned batch", {
    stopifnot(calls == 1L, length(result$pgen) == 4L,
      identical(basename(result$pgen), paste0(c("batch_000003", "batch_000008", "batch_000009", "batch_000010"), ".pgen")))
    for (index in seq_along(result$pgen)) {
      reader <- openPgenReader(result$pgen[[index]], result$pvar[[index]], result$psam[[index]])
      expected <- list(c(1L, 2L, 3L, 5L, 6L, 7L), 8L, c(5L, 6L, 7L), 9L)[[index]]
      stopifnot(identical(reader$sample_ids, ids), identical(reader$variants$position, as.numeric(positions[expected])),
        identical(unname(readPgenWindow(reader, if (index == 4L) "2" else "1", min(positions[expected]), max(positions[expected]))), unname(gt[expected, , drop = FALSE])))
      closePgenReader(reader)
      stopifnot("Retained PVAR rows must preserve exact source metadata" =
        identical(readLines(result$pvar[[index]]), c(pvar_header, pvar_rows[expected])),
        "PSAM must preserve FID, phenotype precision and original bytes" =
        identical(readBin(result$psam[[index]], "raw", n = file.info(result$psam[[index]])$size),
          readBin(opt$psam, "raw", n = file.info(opt$psam)$size)))
    }
  })
  run_named_test("metadata streaming preserves chunk order and rejects changed variant order", {
    expected_paths <- file.path(directory, paste0("expected-", seq_along(worklists), ".pvar"))
    batches <- lapply(worklists, readPgenPreparationWorklist, cis_distance = 10L)
    metadata_with_extra <- file.path(directory, "unselected-extra.pvar")
    writeLines(c(pvar_header, pvar_rows,
      "chrUn_123\t0\textra_zero\tA\tG\t.\tPASS\t.",
      "1\t-1\textra_negative\tA\tG\t.\tPASS\t."), metadata_with_extra)
    counts <- .pgenPreparationWritePvars(metadata_with_extra, batches, expected_paths, chunk_size = 2L)
    stopifnot(identical(counts, c(6, 1, 3, 1)))
    for (index in seq_along(expected_paths)) {
      stopifnot(identical(readLines(result$pvar[[index]]), readLines(expected_paths[[index]])))
      .pgenPreparationValidatePvar(result$pvar[[index]], expected_paths[[index]], chunk_size = 2L)
    }
    altered <- readLines(expected_paths[[1L]])
    altered[c(3L, 4L)] <- altered[c(4L, 3L)]
    bad_order <- file.path(directory, "reordered.pvar")
    writeLines(altered, bad_order)
    expect_error_message(.pgenPreparationValidatePvar(bad_order, expected_paths[[1L]], chunk_size = 2L),
      "variant row order differs")
    expect_error_message(.pgenPreparationValidatePvar(result$pvar[[2L]], expected_paths[[1L]], chunk_size = 2L),
      "variant count differs")
  })
  run_named_test("compressed PSAM is restored as exact plain metadata bytes", {
    bytes <- readBin(opt$psam, "raw", n = file.info(opt$psam)$size)
    compressed <- file.path(directory, "input.psam.gz")
    connection <- gzfile(compressed, "wb")
    writeBin(bytes, connection)
    close(connection)
    restored <- file.path(directory, "restored.psam")
    .pgenPreparationCopyText(compressed, restored)
    stopifnot(identical(readBin(restored, "raw", n = file.info(restored)$size), bytes))
  })
  report_path <- result$reports[[1L]]
  options <- list(pgen = result$pgen[[1L]], pvar = result$pvar[[1L]], psam = result$psam[[1L]],
    worklist = worklists[[1L]], prepared_batch_manifest = report_path, genotype_qc_report = result$qc_report,
    mac = 0, cis_distance = 10L)
  reader <- openPgenReader(options$pgen, options$pvar, options$psam)
  on.exit(closePgenReader(reader), add = TRUE)
  worklist <- readPgenPreparationWorklist(options$worklist, 10L)
  report <- readPgenPreparedBatchReport(options, worklist, reader)
  run_named_test("report binds shared cohort, shard, worklist and windows", {
    stopifnot(report$variant_count == 6L, report$sample_count == 6L,
      identical(report$sample_ids_sha256, digest::digest(ids, algo = "sha256", serializeVersion = 2)),
      identical(report$sample_ids_hash_algorithm, "sha256:R-serialize-v2"),
      !"sample_ids" %in% names(report),
      report$shared_qc$report$variant_count_after == 9L)
    checked <- readPgenBatchQcReport(NULL, report$source_filtered_files, 0, 6L, report$shared_qc$report)
    stopifnot(checked$thresholds$mac == 1L)
  })
  save_report <- function(value) jsonlite::write_json(value, report_path, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  run_named_test("tampered provenance is rejected before fitting", {
    for (field in c("batch_id", "chromosome", "cis_distance", "sample_count", "variant_count")) {
      altered <- report
      altered[[field]] <- if (is.numeric(altered[[field]])) altered[[field]] + 1 else "wrong"
      save_report(altered)
      expect_error_message(readPgenPreparedBatchReport(options, worklist, reader), "Invalid prepared batch report")
    }
    for (role in c("pgen", "pvar", "psam")) {
      altered <- report
      altered$shard_files[[role]]$sha256 <- strrep("0", 64)
      save_report(altered)
      expect_error_message(readPgenPreparedBatchReport(options, worklist, reader), "Invalid prepared batch report")
    }
    altered <- report
    altered$sample_ids_sha256 <- digest::digest(rev(ids), algo = "sha256", serializeVersion = 2)
    save_report(altered)
    expect_error_message(readPgenPreparedBatchReport(options, worklist, reader), "Invalid prepared batch report")
    altered <- report; altered$windows[[1L]]$window_start <- 1L; save_report(altered)
    expect_error_message(readPgenPreparedBatchReport(options, worklist, reader), "Invalid prepared batch report")
    altered <- report; altered$worklist$sha256 <- strrep("0", 64); save_report(altered)
    expect_error_message(readPgenPreparedBatchReport(options, worklist, reader), "Invalid prepared batch report")
    altered <- report; altered$source_filtered_files$pgen$sha256 <- strrep("0", 64); save_report(altered)
    expect_error_message(readPgenPreparedBatchReport(options, worklist, reader), "filtered pgen hash")
    altered <- report; altered$shared_qc$file$sha256 <- strrep("0", 64); save_report(altered)
    expect_error_message(readPgenPreparedBatchReport(options, worklist, reader), "Shared QC report")
    altered <- report; altered$shared_qc$report$sample_count <- 5L; save_report(altered)
    expect_error_message(readPgenPreparedBatchReport(options, worklist, reader), "sample_count")
    for (field in names(report)) {
      altered <- report; altered[[field]] <- NULL; save_report(altered)
      expect_error_message(readPgenPreparedBatchReport(options, worklist, reader), "Invalid prepared batch report")
    }
    without_qc_report <- options
    without_qc_report$genotype_qc_report <- NULL
    for (case in c("original", "changed_threshold", "changed_report_hash")) {
      altered <- report
      if (case == "changed_threshold") altered$shared_qc$report$thresholds$mac <- 999L
      if (case == "changed_report_hash") altered$shared_qc$file$sha256 <- strrep("0", 64)
      save_report(altered)
      expect_error_message(readPgenPreparedBatchReport(without_qc_report, worklist, reader),
        "requires --genotype-qc-report")
    }
    save_report(report)
    options$mac <- 1
    expect_error_message(readPgenPreparedBatchReport(options, worklist, reader), "--mac to 0")
    options$mac <- 0
    options$pgen <- "gs://bucket/unlocalized.pgen"
    expect_error_message(readPgenPreparedBatchReport(options, worklist, reader), "localization error")
  })
  run_named_test("disabled QC uses original cohort with no filter call", {
    opt$mac <- opt$min_call_rate <- 0
    opt$output_dir <- file.path(directory, "disabled")
    disabled <- preparePgenBatches(opt)
    stopifnot(calls == 1L, is.null(disabled$qc_report), identical(disabled$whole_pgen, normalizePath(opt$pgen)),
      is.null(jsonlite::read_json(disabled$reports[[1L]])$shared_qc))
    disabled_options <- options
    disabled_options$pgen <- disabled$pgen[[1L]]
    disabled_options$pvar <- disabled$pvar[[1L]]
    disabled_options$psam <- disabled$psam[[1L]]
    disabled_options$prepared_batch_manifest <- disabled$reports[[1L]]
    disabled_options$genotype_qc_report <- NULL
    disabled_reader <- openPgenReader(disabled_options$pgen, disabled_options$pvar, disabled_options$psam)
    checked <- readPgenPreparedBatchReport(disabled_options, worklist, disabled_reader)
    closePgenReader(disabled_reader)
    stopifnot(is.null(checked$shared_qc))
  })
}
run_preparation_tests()
message("PGEN preparation tests passed")
