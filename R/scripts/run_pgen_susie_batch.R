#!/usr/bin/env Rscript
# Extract a union of cis windows, then run and save each gene separately.
benchmark_library <- Sys.getenv("BENCH_LIBRARY", unset = "")
if (nzchar(benchmark_library)) .libPaths(c(benchmark_library, .libPaths()))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(purrr))
functions_dir <- Sys.getenv("SUSIER_FUNCTIONS_PATH", unset = "/opt/r/lib")
source_files <- vapply(sys.frames(), function(frame) {
  if (is.null(frame$ofile)) "" else frame$ofile
}, character(1))
batch_script_path <- if (any(nzchar(source_files))) tail(source_files[nzchar(source_files)], 1L) else {
  sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][[1L]])
}
pipeline_function_files <- c("ImportFunctions.R", "InitFunctions.R", "SusiePerformance.R",
                   "MemoryHelpers.R", "SusieFunctions.R", "PgenFunctions.R",
                   "SusieOutputFunctions.R", "PgenGeneFunctions.R",
                   "PgenBatchFunctions.R")
for (filename in pipeline_function_files) {
  source(file.path(functions_dir, filename))
}

additionalGenotypesForGene <- function(additional_table, phenotype_ids) {
  if (is.null(additional_table)) return(NULL)
  # Use the old per-gene loader conversion. Keep source row/sample order and
  # missing codes, including the sample columns when there are no feature rows.
  additional_table |>
    dplyr::filter(.data$phenotype_id %in% phenotype_ids) |>
    dplyr::select(-2, -3) |>
    tibble::column_to_rownames("variant_id") |>
    data.matrix()
}

makePgenBatchProvenance <- function(opt, batch_id, batch_files, plink_version) {
  # Cache by local path across roles. Hash scientific inputs and extracted files
  # once per batch, before the gene loop; never hash a chromosome for each gene.
  hashes <- new.env(parent = emptyenv())
  file_record <- function(path) {
    if (is.null(path)) return(list(present = FALSE, sha256 = NA_character_))
    path <- normalizePath(path, mustWork = TRUE)
    if (!exists(path, envir = hashes, inherits = FALSE)) {
      assign(path, digest::digest(file = path, algo = "sha256"), envir = hashes)
    }
    list(present = TRUE, basename = basename(path), size_bytes = file.info(path)$size,
         sha256 = get(path, envir = hashes, inherits = FALSE))
  }
  roles <- c("pgen", "pvar", "psam", "worklist", "phenotype_bed", "qtl_results",
             "covariates", "sample_list", "ancestry", "variant_list", "additional_genotypes")
  inputs <- purrr::map(setNames(roles, roles), ~ file_record(opt[[.x]]))
  scripts_dir <- dirname(normalizePath(batch_script_path, mustWork = TRUE))
  source_paths <- c(file.path(functions_dir, pipeline_function_files),
    file.path(scripts_dir, c("plan_pgen_susie.R", "run_pgen_susie_batch.R",
                            "aggregate_pgen_susie.R")))
  pipeline <- purrr::map(setNames(source_paths, basename(source_paths)), file_record)
  filter_mode <- if (!is.null(opt$variant_list)) "variant_list" else
    if (!is.null(opt$ancestry) && opt$maf > 0) "ancestry_maf" else "none"
  list(batch_id = batch_id, plink_version = plink_version, inputs = inputs,
    extracted_genotypes = purrr::map(setNames(batch_files, c("pgen", "pvar", "psam")), file_record),
    settings = list(maf = opt$maf, variant_filter_mode = filter_mode,
      memory_debug = isTRUE(opt$memory_debug),
      cis_distance = opt$cis_distance, missing_code = opt$missing_code,
      select_top_phenotype_per_cluster = isTRUE(opt$select_top_phenotype_per_cluster),
      top_phenotype_pvalue_column = opt$top_phenotype_pvalue_column,
      base_seed = opt$base_seed,
      seed_algorithm = "susier-gene-outcome-seed-v1:UTF8-byte-length:SHA256-first7hex",
      additional_genotype_scope = "union_of_selected_phenotypes_within_gene",
      plink_threads = 1L, plink_memory_mb = opt$plink_memory_mb),
    source = list(pipeline = pipeline,
      runtime_image = Sys.getenv("SUSIER_RUNTIME_IMAGE", unset = "local"),
      susie_source_sha = Sys.getenv("SUSIER_SUSIER_SOURCE_SHA", unset = "unavailable"),
      eqtlutils_source_sha = Sys.getenv("SUSIER_EQTLUTILS_SOURCE_SHA", unset = "unavailable"),
      r_version = as.character(getRversion())),
    susie_version = as.character(utils::packageVersion("susieR")),
    eqtlutils_version = as.character(utils::packageVersion("eQTLUtils")),
    pgenlibr_version = as.character(utils::packageVersion("pgenlibr")))
}

requireBatchLocalFiles <- function(paths) {
  paths <- unlist(paths, use.names = TRUE)
  if (!length(paths) || anyNA(paths) || any(!nzchar(paths))) {
    stop("Required input paths are missing", call. = FALSE)
  }
  if (any(grepl("^[A-Za-z][A-Za-z0-9+.-]*://", paths))) {
    stop("Input localization error: unresolved cloud URI; supply localized files",
         call. = FALSE)
  }
  if (any(grepl("[\r\n]", paths)) || any(!file.exists(paths)) ||
      any(dir.exists(paths)) || any(file.access(paths, 4L) != 0L)) {
    stop("Batch inputs must be readable local files with single-line paths",
         call. = FALSE)
  }
  invisible(TRUE)
}

decorateBatchTable <- function(table, worklist, gene_id) {
  # molecular_trait_id stays the original outcome ID. Genomic chromosome fields
  # already present in the legacy variant table retain their representation.
  mapping <- worklist |>
    dplyr::filter(.data$gene_id == .env$gene_id) |>
    dplyr::select(batch_id, chromosome, gene_id, outcome_id, phenotype_id)
  table <- table |>
    dplyr::mutate(molecular_trait_id = as.character(molecular_trait_id))
  if ("chromosome" %in% names(table)) {
    mapping <- dplyr::select(mapping, -chromosome)
  }
  mapping <- dplyr::mutate(mapping, molecular_trait_id = outcome_id)
  dplyr::left_join(table, mapping, by = "molecular_trait_id")
}

makePgenBatchExperiment <- function(data) {
  # The upstream constructor drops dimensions for a one-row batch. Keep its
  # sample/feature intersection and ordering, with explicit matrix dimensions.
  samples <- as.data.frame(data$sample_metadata)
  rownames(samples) <- samples$sample_id
  assay <- as.data.frame(dplyr::select(data$expression_matrix, -1, -2, -3))
  rownames(assay) <- assay$phenotype_id
  metadata <- as.data.frame(data$phenotype_meta)
  rownames(metadata) <- metadata$phenotype_id
  shared_samples <- intersect(samples$sample_id, colnames(assay))
  shared_phenotypes <- intersect(metadata$phenotype_id, rownames(assay))
  if (!length(shared_samples) || !length(shared_phenotypes)) {
    stop("No shared phenotype samples or features for this batch", call. = FALSE)
  }
  assay <- as.matrix(assay[shared_phenotypes, shared_samples, drop = FALSE])
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = assay),
    rowData = metadata[shared_phenotypes, , drop = FALSE],
    colData = samples[shared_samples, , drop = FALSE]
  )
  eQTLUtils::subsetSEByColumnValue(se, "qtl_group", "ALL")
}

runPgenSusieBatch <- function(opt) {
  previous_memory_debug <- options(susier.memory_debug = isTRUE(opt$memory_debug))
  on.exit(options(previous_memory_debug), add = TRUE)
  required_names <- c("pgen", "pvar", "psam", "worklist", "phenotype_bed",
                      "qtl_results", "covariates")
  if (any(vapply(opt[required_names], is.null, logical(1)))) {
    stop("Supply all required PGEN, worklist and phenotype input arguments",
         call. = FALSE)
  }
  requireBatchLocalFiles(opt[required_names])
  optional_names <- c("sample_list", "ancestry", "variant_list", "additional_genotypes")
  optional_paths <- opt[optional_names]
  optional_paths <- optional_paths[!vapply(optional_paths, is.null, logical(1))]
  if (length(optional_paths)) requireBatchLocalFiles(optional_paths)
  if (is.null(opt$missing_code) || !opt$missing_code %in% c("NA", "-1")) {
    stop("--missing-code is required and must be NA or -1", call. = FALSE)
  }
  if (is.null(opt$output_dir) || !nzchar(opt$output_dir) ||
      grepl("[\r\n]", opt$output_dir)) {
    stop("--output-dir must be a nonempty local directory", call. = FALSE)
  }
  if (grepl("^[A-Za-z][A-Za-z0-9+.-]*://", opt$output_dir)) {
    stop("Output directory must be local to the task", call. = FALSE)
  }
  valid_integer <- function(value, lower) {
    is.numeric(value) && length(value) == 1L && is.finite(value) &&
      value >= lower && value <= .Machine$integer.max && value == floor(value)
  }
  if (!valid_integer(opt$cis_distance, 0) || !valid_integer(opt$base_seed, 0) ||
      !valid_integer(opt$plink_memory_mb, 64) ||
      length(opt$maf) != 1L || !is.finite(opt$maf) || opt$maf < 0 || opt$maf > 0.5) {
    stop("Invalid cis distance, seed, PLINK memory or MAF", call. = FALSE)
  }
  worklist <- readr::read_tsv(
    opt$worklist, col_types = readr::cols(.default = readr::col_character()),
    na = character(), show_col_types = FALSE
  )
  required_columns <- c("batch_id", "chromosome", "gene_id", "outcome_id",
                         "phenotype_id", "phenotype_pos", "window_start",
                         "window_end", "region")
  if (!all(required_columns %in% names(worklist)) || !nrow(worklist) ||
      anyDuplicated(worklist$outcome_id) || anyDuplicated(worklist$phenotype_id) ||
      anyNA(worklist[required_columns]) ||
      any(!nzchar(as.matrix(worklist[required_columns]))) ||
      dplyr::n_distinct(worklist$batch_id) != 1L ||
      dplyr::n_distinct(worklist$chromosome) != 1L) {
    stop("Invalid batch worklist or duplicate outcome IDs", call. = FALSE)
  }
  numeric_columns <- c("phenotype_pos", "requested_window_start",
                         "requested_window_end", "window_start", "window_end")
  if (!all(numeric_columns %in% names(worklist))) {
    stop("Batch worklist is missing its requested window bounds", call. = FALSE)
  }
  worklist <- worklist |>
    dplyr::mutate(dplyr::across(dplyr::all_of(numeric_columns),
                               ~ suppressWarnings(as.numeric(.x))))
  if (any(!is.finite(as.matrix(worklist[numeric_columns]))) ||
      any(as.matrix(worklist[numeric_columns]) != floor(as.matrix(worklist[numeric_columns]))) ||
      any(worklist$requested_window_start != worklist$phenotype_pos - opt$cis_distance) ||
      any(worklist$requested_window_end != worklist$phenotype_pos + opt$cis_distance) ||
      any(worklist$window_start != pmax(1, worklist$requested_window_start)) ||
      any(worklist$window_end != worklist$requested_window_end)) {
    stop("Invalid batch cis windows for the supplied cis distance", call. = FALSE)
  }
  batch_id <- worklist$batch_id[[1L]]
  message("START batch ", batch_id, ": ", dplyr::n_distinct(worklist$gene_id),
          " genes, ", nrow(worklist), " phenotype fits")
  dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(opt$output_dir, mustWork = TRUE)
  regions_file <- file.path(output_dir, "batch.windows.tsv")
  writePgenBatchRegions(worklist, regions_file)
  pgen_prefix <- file.path(output_dir, "batch.genotypes")
  plink <- opt$plink2
  if (!grepl("/", plink, fixed = TRUE)) plink <- Sys.which(plink)
  if (!nzchar(plink) || !file.exists(plink) || file.access(plink, 1L) != 0L) {
    stop("PLINK 2 executable is not available: ", opt$plink2, call. = FALSE)
  }
  version <- system2(plink, "--version", stdout = TRUE, stderr = TRUE)
  message("PLINK version: ", paste(version, collapse = " "))
  plink_args <- c("--pgen", normalizePath(opt$pgen),
                  "--pvar", normalizePath(opt$pvar),
                  "--psam", normalizePath(opt$psam),
                  "--extract", "range", regions_file, "--make-pgen",
                  "--threads", "1", "--memory", as.character(opt$plink_memory_mb),
                  "--out", pgen_prefix)
  message("EXTRACT batch cis windows with PLINK 2 (1-based inclusive ranges)")
  plink_log <- file.path(output_dir, "plink.stdout.log")
  status <- system2(plink, shQuote(plink_args), stdout = plink_log, stderr = plink_log)
  if (!identical(status, 0L)) {
    detail <- paste(tail(readLines(plink_log, warn = FALSE), 30L), collapse = "\n")
    stop("PLINK multi-region extraction failed (", status, "):\n", detail,
         call. = FALSE)
  }
  batch_files <- paste0(pgen_prefix, c(".pgen", ".pvar", ".psam"))
  requireBatchLocalFiles(batch_files)
  # Keep the complete genotype cohort for ancestry MAF. Only the phenotype
  # metadata uses the inferred sample intersection.
  reader <- openPgenReader(batch_files[[1L]], batch_files[[2L]], batch_files[[3L]])
  on.exit(closePgenReader(reader), add = TRUE)
  message("OPEN batch PGEN: ", reader$raw_sample_count, " samples, ",
          nrow(reader$variants), " variants")
  legacy_options <- list(
    expression_matrix = opt$phenotype_bed, phenotype_list = opt$qtl_results,
    covariates = opt$covariates, genotype_matrix = batch_files[[1L]],
    sample_meta = opt$sample_list, cv_meta = NULL, out_prefix = batch_id,
    cisdistance = opt$cis_distance, n_folds = NULL, VariantList = opt$variant_list,
    MAF = opt$maf, AncestryMetadata = opt$ancestry,
    # Parse the additional table once below, before its ownership columns would
    # be removed by LoadData. Do not allocate a second full batch matrix.
    AdditionalGenotypesBed = NULL,
    select_top_phenotype_per_cluster = isTRUE(opt$select_top_phenotype_per_cluster),
    top_phenotype_pvalue_column = opt$top_phenotype_pvalue_column
  )
  if (is.null(opt$sample_list)) {
    legacy_options$genotype_sample_ids <- reader$sample_ids
  }
  data <- LoadData(legacy_options)
  additional_table <- if (is.null(opt$additional_genotypes)) NULL else
    ImportAdditionalGenotypes(opt$additional_genotypes)
  selected_ids <- as.character(data$phenotype_list$phenotype_id)
  if (!setequal(selected_ids, worklist$phenotype_id)) {
    stop("Batch worklist does not match the phenotypes selected by the current loader",
         call. = FALSE)
  }
  metadata_match <- match(worklist$phenotype_id, data$phenotype_meta$phenotype_id)
  metadata <- data$phenotype_meta[metadata_match, , drop = FALSE]
  actual_positions <- as.numeric(metadata$phenotype_pos)
  requested_starts <- actual_positions - opt$cis_distance
  requested_ends <- actual_positions + opt$cis_distance
  if (anyNA(actual_positions) ||
      any(actual_positions != as.numeric(worklist$phenotype_pos)) ||
      any(normalize_chr(metadata$chromosome) != normalize_chr(worklist$chromosome)) ||
      any(pmax(1, requested_starts) != as.numeric(worklist$window_start)) ||
      any(requested_ends != as.numeric(worklist$window_end))) {
    stop("Batch cis windows differ from the phenotype coordinates or cis distance",
         call. = FALSE)
  }
  se <- makePgenBatchExperiment(data)
  provenance <- makePgenBatchProvenance(opt, batch_id, batch_files, version)
  provenance$settings$sample_selection <- if (is.null(opt$sample_list))
    "genotype_phenotype_intersection" else "provided_list"
  provenance$settings$sample_order <- if (is.null(opt$sample_list)) "psam" else "sample_list"
  provenance$settings$sample_count <- ncol(se)
  provenance$settings$sample_ids_sha256 <- digest::digest(
    as.character(SummarizedExperiment::colData(se)$genotype_id),
    algo = "sha256", serializeVersion = 2
  )
  compact <- purrr::map(unique(worklist$gene_id), function(gene_id) {
    selected <- dplyr::filter(worklist, .data$gene_id == .env$gene_id)
    additional_genotypes <- additionalGenotypesForGene(additional_table, selected$phenotype_id)
    outcomes <- purrr::map(seq_len(nrow(selected)), function(i) {
      list(phenotype_id = selected$phenotype_id[[i]], se = se,
           covariates = data$covariates_matrix, cis_distance = opt$cis_distance,
           ancestry_data = data$AncestryDf, MAF = data$MAF_threshold,
           variant_list = data$variant_list,
           additional_genotypes = additional_genotypes,
           region = selected$region[[i]])
    })
    names(outcomes) <- selected$outcome_id
    result <- runPgenGene(gene_id, outcomes, reader, output_dir, opt$base_seed,
                           opt$missing_code, provenance)
    result$summary <- decorateBatchTable(result$summary, worklist, gene_id)
    result$credible_sets <- decorateBatchTable(result$credible_sets, worklist, gene_id)
    result$diagnostics <- decorateBatchTable(result$diagnostics, worklist, gene_id)
    result$index <- dplyr::mutate(
      result$index, batch_id = batch_id, chromosome = worklist$chromosome[[1L]]
    )
    message("COMPLETE gene ", gene_id)
    result
  })
  # Old empty credible-set schemas have numeric ID columns. Empty rows must not
  # coerce populated tables or cause an otherwise valid mixed batch to fail.
  bind_compact <- function(name) {
    tables <- purrr::map(compact, name)
    nonempty <- purrr::keep(tables, ~ nrow(.x) > 0L)
    if (length(nonempty)) dplyr::bind_rows(nonempty) else tables[[1L]]
  }
  arrow::write_parquet(bind_compact("summary"),
                        file.path(output_dir, "batch.susie.parquet"))
  arrow::write_parquet(bind_compact("credible_sets"),
                        file.path(output_dir, "batch.credible_sets.parquet"))
  readr::write_tsv(bind_compact("index"),
                    file.path(output_dir, "batch.outcome_index.tsv"))
  readr::write_tsv(bind_compact("diagnostics"),
                    file.path(output_dir, "batch.covariate_explained_variants.tsv"))
  message("COMPLETE batch ", batch_id)
  invisible(TRUE)
}

batchOption <- function(name, type = "character", default = NULL, help = name) {
  optparse::make_option(paste0("--", name), type = type, default = default,
                        help = help, dest = gsub("-", "_", name, fixed = TRUE))
}
batch_options <- list(
  batchOption("pgen"), batchOption("pvar"), batchOption("psam"),
  batchOption("worklist"), batchOption("phenotype-bed"), batchOption("qtl-results"),
  batchOption("covariates"), batchOption("sample-list"), batchOption("missing-code"),
  batchOption("output-dir"), batchOption("cis-distance", "integer", 1000000L),
  batchOption("base-seed", "integer", 1L), batchOption("maf", "double", 0),
  batchOption("ancestry"), batchOption("variant-list"),
  batchOption("additional-genotypes"), batchOption("plink2", default = "plink2"),
  batchOption("plink-memory-mb", "integer", 1024L),
  optparse::make_option("--memory-debug", action = "store_true", default = FALSE,
                        dest = "memory_debug", help = "Print detailed memory diagnostics"),
  optparse::make_option("--select-top-phenotype-per-cluster", action = "store_true",
                        default = FALSE, dest = "select_top_phenotype_per_cluster"),
  batchOption("top-phenotype-pvalue-column", default = "qval")
)
if (sys.nframe() == 0L) {
  parser <- optparse::OptionParser(option_list = batch_options)
  runPgenSusieBatch(optparse::parse_args(parser))
}
