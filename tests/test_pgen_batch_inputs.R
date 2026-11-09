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

check_mac_cli_and_validation <- function() {
  parser <- optparse::OptionParser(option_list = batch_options)
  defaults <- optparse::parse_args(parser, args = character(), print_help_and_exit = FALSE)
  stopifnot("The MAC filter must be off by default" = identical(defaults$mac, 0))
  enabled <- optparse::parse_args(parser, args = c("--mac", "5"), print_help_and_exit = FALSE)
  fractional <- optparse::parse_args(parser, args = c("--mac", "1.5"), print_help_and_exit = FALSE)
  stopifnot(identical(enabled$mac, 5), identical(fractional$mac, 1.5))

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
if (Sys.getenv("PGEN_INPUT_CHECK", "all") %in% c("all", "provenance")) check_batch_provenance()
message("PGEN batch input tests passed")
