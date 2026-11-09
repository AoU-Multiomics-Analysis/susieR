# Input ownership and provenance checks without PLINK or model fitting.
benchmark_library <- Sys.getenv("BENCH_LIBRARY", unset = "")
if (nzchar(benchmark_library)) .libPaths(c(benchmark_library, .libPaths()))
repo_dir <- normalizePath(getwd())
scripts_dir <- Sys.getenv("SUSIER_SCRIPTS_PATH", file.path(repo_dir, "R", "scripts"))
source(file.path(scripts_dir, "run_pgen_susie_batch.R"))

check_additional_ownership <- function() {
  # Catch loss of ownership, row/sample reordering, sentinel changes, and NULL
  # replacement of a zero-row matrix (which would skip sample intersection).
  raw <- tibble::tibble(variant_id = c("b2", "a", "b1"),
    phenotype_id = c("outB2", "outA", "outB1"), feature_type = "SV",
    `109` = c(NA_real_, 2, -1), `102` = c(1, -1, 0))
  b <- additionalGenotypesForGene(raw, c("outB1", "outB2"))
  stopifnot(identical(rownames(b), c("b2", "b1")),
    identical(colnames(b), c("109", "102")),
    identical(unname(b), matrix(c(NA_real_, -1, 1, 0), nrow = 2L)))
  empty <- additionalGenotypesForGene(raw, "outC")
  stopifnot(is.matrix(empty), identical(dim(empty), c(0L, 2L)),
    identical(colnames(empty), c("109", "102")),
    is.null(additionalGenotypesForGene(NULL, "outC")))
  gt <- matrix(0, nrow = 1L, ncol = 3L,
    dimnames = list("original", c("102", "108", "109")))
  stopifnot(identical(colnames(MergeAdditionalGenotypes(gt, empty)), c("102", "109")))
}

check_batch_provenance <- function() {
  # Catch omitted optional files, role collisions, settings, and source code.
  directory <- tempfile("batch-provenance-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  roles <- c("pgen", "pvar", "psam", "worklist", "phenotype_bed", "qtl_results",
             "covariates", "sample_list", "additional_genotypes")
  opt <- setNames(lapply(roles, function(role) {
    path <- file.path(directory, role, "same-name")
    dir.create(dirname(path))
    writeLines(role, path)
    path
  }), roles)
  opt <- c(opt, list(maf = 0, cis_distance = 100L, missing_code = "NA",
    base_seed = 82L, select_top_phenotype_per_cluster = FALSE,
    top_phenotype_pvalue_column = "qval", plink_memory_mb = 1024L))
  batch_files <- unlist(opt[c("pgen", "pvar", "psam")], use.names = FALSE)
  a <- makePgenBatchProvenance(opt, "batch1", batch_files, "PLINK test")
  stopifnot(identical(a$settings$memory_debug, FALSE),
    identical(a$settings$mac, 0),
    identical(a$settings$min_call_rate, 0),
    identical(a$settings$sample_scope, "disabled"),
    identical(a$settings$genotype_qc_count_source, "disabled"),
    identical(a$settings$mac_sample_scope, "all_input_pgen_samples"),
    identical(a$settings$mac_include_nonfounders, TRUE))
  opt$memory_debug <- TRUE
  debug <- makePgenBatchProvenance(opt, "batch1", batch_files, "PLINK test")
  stopifnot(identical(debug$settings$memory_debug, TRUE))
  manifest <- file.path(directory, "provenance.json")
  jsonlite::write_json(a, manifest, auto_unbox = TRUE, na = "null")
  serialized <- jsonlite::read_json(manifest)
  stopifnot("Absent input hashes must be JSON null" = is.null(serialized$inputs$ancestry$sha256))
  stopifnot(all(roles %in% names(a$inputs)),
    a$inputs$additional_genotypes$present, !a$inputs$ancestry$present,
    !a$inputs$variant_list$present,
    identical(a$inputs$pgen$sha256, digest::digest(file = opt$pgen, algo = "sha256")),
    !identical(a$inputs$pgen$sha256, a$inputs$pvar$sha256),
    identical(a$settings$maf, 0), identical(a$settings$cis_distance, 100L),
    identical(a$settings$missing_code, "NA"), identical(a$settings$base_seed, 82L),
    !a$settings$select_top_phenotype_per_cluster,
    identical(a$settings$top_phenotype_pvalue_column, "qval"),
    identical(a$source$pipeline$run_pgen_susie_batch.R$sha256,
      digest::digest(file = file.path(scripts_dir, "run_pgen_susie_batch.R"), algo = "sha256")),
    nzchar(a$source$pipeline$ImportFunctions.R$sha256),
    identical(a$source$runtime_image, Sys.getenv("SUSIER_RUNTIME_IMAGE", "local")))
  writeLines("changed additional genotypes", opt$additional_genotypes)
  opt$maf <- .1
  b <- makePgenBatchProvenance(opt, "batch1", batch_files, "PLINK test")
  stopifnot(!identical(a$inputs$additional_genotypes$sha256,
                      b$inputs$additional_genotypes$sha256),
    identical(b$settings$maf, .1), identical(a$inputs$covariates, b$inputs$covariates),
    identical(b$settings$variant_filter_mode, "none"))
  opt$ancestry <- opt$covariates
  ancestry <- makePgenBatchProvenance(opt, "batch1", batch_files, "PLINK test")
  stopifnot(ancestry$inputs$ancestry$present,
    identical(ancestry$settings$variant_filter_mode, "ancestry_maf"))
  opt$variant_list <- opt$sample_list
  variants <- makePgenBatchProvenance(opt, "batch1", batch_files, "PLINK test")
  stopifnot(variants$inputs$variant_list$present,
    identical(variants$settings$variant_filter_mode, "variant_list"))
  opt$mac <- 5
  mac_variants <- makePgenBatchProvenance(opt, "batch1", batch_files, "PLINK test")
  stopifnot(identical(mac_variants$settings$mac, 5),
    identical(mac_variants$settings$variant_filter_mode, "mac+variant_list"))
  opt$variant_list <- NULL
  mac_ancestry <- makePgenBatchProvenance(opt, "batch1", batch_files, "PLINK test")
  stopifnot(identical(mac_ancestry$settings$variant_filter_mode, "mac+ancestry_maf"))
  opt$ancestry <- NULL
  mac_only <- makePgenBatchProvenance(opt, "batch1", batch_files, "PLINK test")
  stopifnot(identical(mac_only$settings$variant_filter_mode, "mac"))
}

check_shared_qc_provenance <- function() {
  # A report must identify these exact prepared files before a fit can use it.
  directory <- tempfile("batch-shared-qc-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  roles <- c("pgen", "pvar", "psam", "worklist", "phenotype_bed", "qtl_results",
             "covariates")
  opt <- setNames(lapply(roles, function(role) {
    path <- file.path(directory, role)
    writeLines(role, path)
    path
  }), roles)
  opt <- c(opt, list(maf = 0, mac = 0, cis_distance = 100L, missing_code = "NA",
    base_seed = 82L, select_top_phenotype_per_cluster = FALSE,
    top_phenotype_pvalue_column = "qval", plink_memory_mb = 1024L))
  batch_files <- unlist(opt[c("pgen", "pvar", "psam")], use.names = FALSE)
  baseline <- makePgenBatchProvenance(opt, "batch1", batch_files, "PLINK test")
  counts <- file.path(directory, "counts.tsv.gz")
  writeLines("QC counts", counts)
  record <- function(path) list(present = TRUE, basename = basename(path),
    size_bytes = file.info(path)$size, sha256 = digest::digest(file = path, algo = "sha256"))
  source_records <- purrr::map(setNames(c("pgen", "pvar", "psam"), c("pgen", "pvar", "psam")),
    function(role) {
      path <- file.path(directory, paste0("original.", role))
      writeLines(paste("Original genotype", role), path)
      record(path)
    })
  report <- list(schema_version = 1L, count_source = "pgen_hardcalls",
    sample_scope = "all_input_pgen_samples", include_nonfounders = TRUE,
    sample_count = 10L, variant_count_before = 12L, variant_count_after = 8L,
    thresholds = list(mac = 5L, min_call_rate = .95),
    source_files = source_records,
    filtered_files = baseline$inputs[c("pgen", "pvar", "psam")],
    counts_file = record(counts), plink_version = "PLINK test")
  opt$genotype_qc_report <- file.path(directory, "genotype_qc.json")
  save_report <- function(value) jsonlite::write_json(value, opt$genotype_qc_report,
    auto_unbox = TRUE, na = "null", pretty = TRUE)
  save_report(report)
  actual <- makePgenBatchProvenance(opt, "batch1", batch_files, "PLINK test")
  stopifnot(identical(actual$genotype_qc$thresholds$mac, 5L),
    identical(actual$genotype_qc$source_files$pgen$sha256, source_records$pgen$sha256),
    !identical(actual$genotype_qc$source_files$pgen$sha256, actual$inputs$pgen$sha256),
    identical(actual$inputs$pgen$sha256, baseline$inputs$pgen$sha256),
    identical(actual$settings$mac, 5L), identical(actual$settings$min_call_rate, .95),
    identical(actual$settings$variant_filter_mode, "mac+call_rate"),
    identical(actual$settings$sample_scope, "all_input_pgen_samples"),
    identical(actual$settings$genotype_qc_count_source, "pgen_hardcalls"),
    isTRUE(actual$settings$include_nonfounders),
    actual$inputs$genotype_qc_report$present,
    identical(actual$inputs$genotype_qc_report$sha256,
      digest::digest(file = opt$genotype_qc_report, algo = "sha256")),
    nzchar(actual$source$pipeline$PgenQcFunctions.R$sha256),
    nzchar(actual$source$pipeline$filter_pgen_susie.R$sha256))
  opt$maf <- .1
  opt$ancestry <- opt$covariates
  ancestry <- makePgenBatchProvenance(opt, "batch1", batch_files, "PLINK test")
  stopifnot(identical(ancestry$settings$variant_filter_mode, "mac+call_rate+ancestry_maf"))
  opt$variant_list <- opt$covariates
  variants <- makePgenBatchProvenance(opt, "batch1", batch_files, "PLINK test")
  stopifnot(identical(variants$settings$variant_filter_mode, "mac+call_rate+variant_list"))
  opt$variant_list <- opt$ancestry <- NULL
  report$thresholds$mac <- 0L
  save_report(report)
  call_rate_only <- makePgenBatchProvenance(opt, "batch1", batch_files, "PLINK test")
  stopifnot(identical(call_rate_only$settings$variant_filter_mode, "call_rate"))
  report$thresholds <- list(mac = 5L, min_call_rate = 0)
  save_report(report)
  mac_only <- makePgenBatchProvenance(opt, "batch1", batch_files, "PLINK test")
  stopifnot(identical(mac_only$settings$variant_filter_mode, "mac"))
  report$thresholds <- list(mac = 0L, min_call_rate = 0)
  report$variant_count_before <- 8L
  save_report(report)
  disabled <- makePgenBatchProvenance(opt, "batch1", batch_files, "PLINK test")
  stopifnot(identical(disabled$settings$variant_filter_mode, "none"))
  report$thresholds <- list(mac = 5L, min_call_rate = .95)
  report$variant_count_before <- 12L
  save_report(report)
  matched_cohort <- makePgenBatchProvenance(opt, "batch1", batch_files, "PLINK test", 10L)
  stopifnot(identical(matched_cohort$genotype_qc$sample_count, 10L))
  wrong_cohort <- tryCatch(
    makePgenBatchProvenance(opt, "batch1", batch_files, "PLINK test", 9L), error = identity)
  stopifnot(inherits(wrong_cohort, "error"),
    grepl("sample_count", conditionMessage(wrong_cohort)))

  expect_invalid <- function(value, pattern = "genotype QC report") {
    save_report(value)
    error <- tryCatch(makePgenBatchProvenance(opt, "batch1", batch_files, "PLINK test"),
      error = identity)
    stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error), ignore.case = TRUE))
  }
  # Removing validation, accepting a stale report, or accepting invalid cohort
  # counts must fail these cases even when the report parses as valid JSON.
  for (role in c("pgen", "pvar", "psam")) {
    wrong <- report
    wrong$filtered_files[[role]]$sha256 <- paste(rep("0", 64), collapse = "")
    expect_invalid(wrong, paste0("filtered.*", role, ".*hash|", role, ".*sha256"))
  }
  for (field in names(report)) {
    missing <- report
    missing[[field]] <- NULL
    expect_invalid(missing)
  }
  invalid_fields <- list(schema_version = 2L, count_source = "pvar_info",
    sample_scope = "phenotype_intersection", include_nonfounders = FALSE,
    sample_count = 0L, variant_count_before = 7L, variant_count_after = 13L,
    plink_version = "")
  for (field in names(invalid_fields)) {
    invalid <- report
    invalid[[field]] <- invalid_fields[[field]]
    expect_invalid(invalid)
  }
  for (field in c("sample_count", "variant_count_before", "variant_count_after")) {
    invalid <- report
    invalid[[field]] <- 1.5
    expect_invalid(invalid)
  }
  for (thresholds in list(list(mac = -1, min_call_rate = .95),
    list(mac = 1.5, min_call_rate = .95), list(mac = 0, min_call_rate = 1.1),
    list(mac = 0, min_call_rate = NA_real_), list(mac = Inf, min_call_rate = .95),
    list(mac = "5", min_call_rate = .95))) {
    invalid <- report
    invalid$thresholds <- thresholds
    expect_invalid(invalid)
  }
  save_report(report)
  writeLines(sub('"mac": 5', '"mac": 1e999',
    readLines(opt$genotype_qc_report), fixed = TRUE), opt$genotype_qc_report)
  nonfinite <- tryCatch(makePgenBatchProvenance(opt, "batch1", batch_files, "PLINK test"),
    error = identity)
  stopifnot(inherits(nonfinite, "error"),
    grepl("genotype QC report", conditionMessage(nonfinite)))
  for (files in c("source_files", "filtered_files")) {
    for (field in c("present", "basename", "size_bytes", "sha256")) {
      invalid <- report
      invalid[[files]]$pgen[[field]] <- NULL
      expect_invalid(invalid)
    }
  }
  invalid <- report
  invalid$counts_file$present <- FALSE
  expect_invalid(invalid)
  invalid <- report
  invalid$thresholds <- list(mac = 0L, min_call_rate = 0)
  expect_invalid(invalid)
  invalid <- report
  invalid$filtered_files$pgen$size_bytes <- 99L
  expect_invalid(invalid)
  writeLines("{not json", opt$genotype_qc_report)
  error <- tryCatch(makePgenBatchProvenance(opt, "batch1", batch_files, "PLINK test"),
    error = identity)
  stopifnot(inherits(error, "error"), grepl("genotype QC report", conditionMessage(error)))
  save_report(report)
  opt$mac <- 1
  expect_invalid(report, "--mac.*0|repeat.*filter")
  opt$mac <- 0
  opt$genotype_qc_report <- NULL
  omitted <- makePgenBatchProvenance(opt, "batch1", batch_files, "PLINK test")
  stopifnot(is.null(omitted$genotype_qc), identical(omitted$settings$mac, 0),
    identical(omitted$settings$min_call_rate, 0),
    !omitted$inputs$genotype_qc_report$present)
}

check_mac_cli_and_validation <- function() {
  parser <- optparse::OptionParser(option_list = batch_options)
  defaults <- optparse::parse_args(parser, args = character(), print_help_and_exit = FALSE)
  stopifnot("The MAC filter must be off by default" = identical(defaults$mac, 0))
  enabled <- optparse::parse_args(parser, args = c("--mac", "5"), print_help_and_exit = FALSE)
  fractional <- optparse::parse_args(parser, args = c("--mac", "1.5"), print_help_and_exit = FALSE)
  stopifnot(identical(enabled$mac, 5), identical(fractional$mac, 1.5))
  qc <- optparse::parse_args(parser, args = c("--genotype-qc-report", "qc.json"),
    print_help_and_exit = FALSE)
  stopifnot(is.null(defaults$genotype_qc_report), identical(qc$genotype_qc_report, "qc.json"))

  directory <- tempfile("batch-mac-inputs-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  roles <- c("pgen", "pvar", "psam", "worklist", "phenotype_bed", "qtl_results", "covariates")
  opt <- setNames(lapply(roles, function(role) {
    path <- file.path(directory, role)
    writeLines("fixture", path)
    path
  }), roles)
  opt <- c(opt, list(missing_code = "NA", output_dir = file.path(directory, "output"),
    cis_distance = 10L, base_seed = 1L, plink_memory_mb = 1024L, maf = 0))
  readr::write_tsv(tibble::tibble(batch_id = "batch1", chromosome = "1", gene_id = "gene1",
    outcome_id = "p1", phenotype_id = "p1", phenotype_pos = 100L,
    requested_window_start = 90L, requested_window_end = 110L,
    window_start = 90L, window_end = 110L, region = "1:90-110"), opt$worklist)
  capture <- file.path(directory, "plink.args")
  opt$plink2 <- file.path(directory, "record-plink")
  # Replace only the external executable. The real runner builds and quotes
  # the command; the deliberate exit stops before opening dummy PGEN files.
  writeLines(c("#!/usr/bin/env bash", "set -euo pipefail",
    'if [[ "${1-}" == --version ]]; then echo "PLINK recording fixture"; exit 0; fi',
    paste0("printf '%s\\n' \"$@\" > ", shQuote(capture)), "exit 1"), opt$plink2)
  Sys.chmod(opt$plink2, mode = "0755")
  failure <- function(options) tryCatch({
    runPgenSusieBatch(options)
    stop("Expected the batch to stop", call. = FALSE)
  }, error = identity)
  report <- opt
  report$genotype_qc_report <- "gs://bucket/qc.json"
  error <- failure(report)
  stopifnot(grepl("localization error", conditionMessage(error)), !file.exists(capture))
  report$genotype_qc_report <- opt$pgen
  report$mac <- 1
  error <- failure(report)
  stopifnot(grepl("--mac to 0", conditionMessage(error)), !file.exists(capture),
    !dir.exists(opt$output_dir))
  for (value in list(-1, 1.5, NA_real_, NaN, Inf, "5", c(1, 2), .Machine$integer.max + 1)) {
    invalid <- opt
    invalid$mac <- value
    error <- failure(invalid)
    stopifnot(inherits(error, "error"), grepl("MAC", conditionMessage(error)),
      !file.exists(capture), !dir.exists(opt$output_dir))
  }
  for (value in list(NULL, 0, 5)) {
    current <- opt
    current$mac <- value
    error <- failure(current)
    stopifnot(grepl("PLINK multi-region extraction failed", conditionMessage(error)),
      file.exists(capture))
    arguments <- readLines(capture)
    if (is.null(value) || value == 0) {
      stopifnot(!"--mac" %in% arguments, !"--nonfounders" %in% arguments)
    } else {
      index <- match("--mac", arguments)
      stopifnot(!is.na(index), identical(arguments[[index + 1L]], "5"),
        sum(arguments == "--mac") == 1L, "--nonfounders" %in% arguments)
    }
    unlink(capture)
    unlink(opt$output_dir, recursive = TRUE)
  }
}
check_memory_debug_cli <- function() {
  parser <- optparse::OptionParser(option_list = batch_options)
  defaults <- optparse::parse_args(parser, args = character(), print_help_and_exit = FALSE)
  stopifnot("Memory diagnostics must be off by default" =
              identical(defaults$memory_debug, FALSE))
  enabled <- optparse::parse_args(parser, args = "--memory-debug", print_help_and_exit = FALSE)
  stopifnot(identical(enabled$memory_debug, TRUE))
  original <- options("susier.memory_debug")
  on.exit(options(original), add = TRUE)
  for (flag in c(FALSE, TRUE)) {
    options(susier.memory_debug = !flag)
    observed <- NULL
    failed <- tryCatch(withCallingHandlers(
      runPgenSusieBatch(list(memory_debug = flag)),
      error = function(condition) observed <<- getOption("susier.memory_debug")
    ), error = function(condition) condition)
    stopifnot(inherits(failed, "error"), identical(observed, flag),
              identical(getOption("susier.memory_debug"), !flag))
  }
}
check_memory_debug_cli()
check_mac_cli_and_validation()
if (Sys.getenv("PGEN_INPUT_CHECK", "all") %in% c("all", "additional")) check_additional_ownership()
if (Sys.getenv("PGEN_INPUT_CHECK", "all") %in% c("all", "provenance")) {
  check_batch_provenance()
  check_shared_qc_provenance()
}
message("PGEN batch input tests passed")
