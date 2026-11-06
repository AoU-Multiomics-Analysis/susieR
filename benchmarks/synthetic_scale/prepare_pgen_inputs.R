# Create matched chr1 inputs for the PGEN loading benchmark.
# Usage: BENCH_LIBRARY=/path/to/lib Rscript prepare_pgen_inputs.R INPUT_DIR OUTPUT_DIR
prepare_pgen_inputs <- function(input_dir, output_dir) {
  required_packages <- c("eQTLUtils", "pgenlibr", "Rsamtools", "dplyr", "readr",
                         "tibble", "stringr", "jsonlite")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages) > 0L) {
    stop("Missing preparation packages: ", paste(missing_packages, collapse = ", "),
         call. = FALSE)
  }
  input_dir <- normalizePath(input_dir, mustWork = TRUE)
  copy_names <- c("phenotypes.bed", "covariates.tsv", "samples.tsv",
                  "permutations.tsv", "causal_variants.tsv", "fixture.json")
  source_paths <- file.path(input_dir, c("dosages.tsv.bgz", "dosages.tsv.bgz.tbi",
                                        copy_names))
  if (any(!file.exists(source_paths)) || any(file.access(source_paths, 4L) != 0L)) {
    stop("The source fixture has missing or unreadable files", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  if (identical(input_dir, output_dir)) {
    stop("The output directory must differ from the source directory", call. = FALSE)
  }
  output_names <- c("dosages.tsv", "dosages.tsv.bgz", "dosages.tsv.bgz.tbi",
                    "synthetic.pgen", "synthetic.pvar", "synthetic.psam",
                    "missing_normalization_masks.tsv", "pgen_preparation.json",
                    copy_names)
  if (any(file.exists(file.path(output_dir, output_names)))) {
    stop("Preparation output files already exist; use a new output directory",
         call. = FALSE)
  }

  message("Load the source chr1 region with eQTLUtils")
  genotype <- eQTLUtils::extractGenotypeMatrixFromDosage(
    chr = "1", start = 1L, end = 536870912L,
    dosage_file = file.path(input_dir, "dosages.tsv.bgz")
  )
  fixture <- jsonlite::read_json(file.path(input_dir, "fixture.json"),
                                 simplifyVector = TRUE)
  if (!is.matrix(genotype) || nrow(genotype) == 0L || ncol(genotype) == 0L ||
      nrow(genotype) != fixture$variants || ncol(genotype) != fixture$samples) {
    stop("The loaded chr1 matrix does not match the source fixture dimensions",
         call. = FALSE)
  }
  sample_ids <- colnames(genotype)
  variant_ids <- rownames(genotype)
  if (is.null(sample_ids) || anyNA(sample_ids) || any(!nzchar(sample_ids)) ||
      anyDuplicated(sample_ids) > 0L || any(grepl("[[:space:]]", sample_ids))) {
    stop("The source matrix requires unique sample IDs without whitespace",
         call. = FALSE)
  }
  if (is.null(variant_ids) || anyNA(variant_ids) || anyDuplicated(variant_ids) > 0L) {
    stop("The source matrix requires unique variant IDs", call. = FALSE)
  }
  parts <- stringr::str_split_fixed(variant_ids, "_", 4L)
  variant_metadata <- tibble::tibble(
    chromosome = stringr::str_remove(parts[, 1L], "^chr"),
    position = parts[, 2L],
    variant_id = variant_ids,
    ref = parts[, 3L],
    alt = parts[, 4L]
  )
  positions <- suppressWarnings(as.numeric(variant_metadata$position))
  if (any(variant_metadata$chromosome != "1") || anyNA(positions) ||
      any(!is.finite(positions)) || any(positions < 1) ||
      any(positions != floor(positions)) || is.unsorted(positions) ||
      any(!nzchar(variant_metadata$ref)) || any(!nzchar(variant_metadata$alt)) ||
      any(grepl(",", variant_metadata$alt, fixed = TRUE))) {
    stop("The source must contain sorted chr1 biallelic variants", call. = FALSE)
  }
  n_samples <- ncol(genotype)
  n_variants <- nrow(genotype)
  legacy_nonzero <- normalized_nonzero <- logical(n_variants)
  original_na <- original_sentinel <- numeric(n_variants)
  observed_values_equal <- missing_masks_equal <- TRUE
  hardcall_counts <- c(`0` = 0, `1` = 0, `2` = 0, missing = 0)

  tsv_path <- file.path(output_dir, "dosages.tsv")
  tsv_connection <- pvar_connection <- writer <- pgen_reader <- pvar_reader <- NULL
  on.exit({
    if (!is.null(pgen_reader)) pgenlibr::ClosePgen(pgen_reader)
    if (!is.null(pvar_reader)) pgenlibr::ClosePvar(pvar_reader)
    if (!is.null(writer)) try(pgenlibr::ClosePgenWriter(writer), silent = TRUE)
    if (!is.null(tsv_connection)) close(tsv_connection)
    if (!is.null(pvar_connection)) close(pvar_connection)
    # The uncompressed table is temporary and can be large.
    if (file.exists(tsv_path)) unlink(tsv_path)
  }, add = TRUE)
  tsv_connection <- file(tsv_path, open = "wt")
  pvar_connection <- file(file.path(output_dir, "synthetic.pvar"), open = "wt")
  writeLines(paste(c("CHROM", "POS", "REF", "ALT", sample_ids), collapse = "\t"),
             tsv_connection)
  writeLines("#CHROM\tPOS\tID\tREF\tALT", pvar_connection)
  readr::write_tsv(tibble::tibble(`#FID` = "0", IID = sample_ids),
                   file.path(output_dir, "synthetic.psam"))
  writer <- pgenlibr::NewPgenWriter(
    file.path(output_dir, "synthetic.pgen"), n_samples, n_variants + 4L,
    nonref_flags = FALSE
  )
  append_pgen_variant <- function(chromosome, position, variant_id, ref, alt, values) {
    pgenlibr::AppendBiallelic(writer, values)
    writeLines(paste(chromosome, position, variant_id, ref, alt, sep = "\t"),
               pvar_connection)
  }

  # Only PGEN has these flanks. Tabix requires one block per chromosome.
  append_pgen_variant("2", "1", "chr2_1_A_G", "A", "G",
                       rep(c(0, 1, 2), length.out = n_samples))
  append_pgen_variant("2", "2", "chr2_2_A_G", "A", "G",
                       rep(c(2, 0, 1, NA_real_), length.out = n_samples))
  message("Write normalized TSV and PGEN chr1 hardcalls")
  for (i in seq_len(n_variants)) {
    # Normalize one row at a time to avoid another full-size matrix or mask.
    original <- genotype[i, ]
    legacy_nonzero[[i]] <- sum(round(original, 0), na.rm = TRUE) != 0
    original_na[[i]] <- sum(is.na(original))
    original_sentinel[[i]] <- sum(original == -1, na.rm = TRUE)
    original_missing <- is.na(original) | original == -1
    values <- original
    values[values == -1] <- NA_real_
    if (any(!is.na(values) & !values %in% c(0, 1, 2))) {
      stop("Invalid hardcall at variant ", variant_ids[[i]], call. = FALSE)
    }
    normalized_nonzero[[i]] <- sum(round(values, 0), na.rm = TRUE) != 0
    observed_values_equal <- observed_values_equal &&
      identical(original[!original_missing], values[!original_missing])
    missing_masks_equal <- missing_masks_equal &&
      identical(original_missing, is.na(values))
    hardcall_counts <- hardcall_counts + c(
      sum(values == 0, na.rm = TRUE), sum(values == 1, na.rm = TRUE),
      sum(values == 2, na.rm = TRUE), sum(is.na(values))
    )
    meta <- variant_metadata[i, ]
    append_pgen_variant(meta$chromosome, meta$position, meta$variant_id,
                         meta$ref, meta$alt, values)
    dosage_fields <- as.character(values)
    dosage_fields[is.na(values)] <- "NA"
    writeLines(paste(c(meta$chromosome, meta$position, meta$ref, meta$alt,
                       dosage_fields), collapse = "\t"), tsv_connection)
    if (i == 1L || i %% 1000L == 0L || i == n_variants) {
      message("Prepared ", i, "/", n_variants, " chr1 variants")
    }
  }
  append_pgen_variant("2", "3", "chr2_3_A_G", "A", "G",
                       rep(c(1, 2, 0), length.out = n_samples))
  append_pgen_variant("2", "4", "chr2_4_A_G", "A", "G",
                       rep(c(NA_real_, 2, 1, 0), length.out = n_samples))
  pgenlibr::ClosePgenWriter(writer)
  writer <- NULL
  close(tsv_connection)
  tsv_connection <- NULL
  close(pvar_connection)
  pvar_connection <- NULL
  message("Check every chr1 PGEN hardcall in 100-variant blocks")
  pvar_reader <- pgenlibr::NewPvar(file.path(output_dir, "synthetic.pvar"))
  pgen_reader <- pgenlibr::NewPgen(file.path(output_dir, "synthetic.pgen"),
                                  pvar = pvar_reader)
  if (pgenlibr::GetRawSampleCt(pgen_reader) != n_samples ||
      pgenlibr::GetVariantCt(pgen_reader) != n_variants + 4L) {
    stop("PGEN roundtrip dimensions do not match", call. = FALSE)
  }
  roundtrip_variants_checked <- 0L
  roundtrip_cells_checked <- 0
  for (first_row in seq.int(1L, n_variants, by = 100L)) {
    rows <- seq.int(first_row, min(first_row + 99L, n_variants))
    expected <- genotype[rows, , drop = FALSE]
    expected[expected == -1] <- NA_real_
    actual <- t(pgenlibr::ReadIntList(pgen_reader, as.integer(rows + 2L)))
    storage.mode(expected) <- storage.mode(actual) <- "double"
    dimnames(expected) <- dimnames(actual) <- NULL
    if (!identical(actual, expected)) {
      stop("PGEN hardcall roundtrip failed for chr1 rows ", first_row, "-",
           tail(rows, 1L), call. = FALSE)
    }
    roundtrip_variants_checked <- roundtrip_variants_checked + length(rows)
    roundtrip_cells_checked <- roundtrip_cells_checked + length(expected)
    if (roundtrip_variants_checked %% 5000L == 0L ||
        roundtrip_variants_checked == n_variants) {
      message("Checked ", roundtrip_variants_checked, "/", n_variants,
              " chr1 variants exactly")
    }
  }
  pgenlibr::ClosePgen(pgen_reader)
  pgen_reader <- NULL
  pgenlibr::ClosePvar(pvar_reader)
  pvar_reader <- NULL
  rm(expected, actual)
  rm(genotype)
  invisible(gc())

  if (!observed_values_equal || !missing_masks_equal) {
    stop("Missing-value normalization changed observed data", call. = FALSE)
  }
  normalization_masks <- variant_metadata |>
    dplyr::transmute(
      variant_id,
      legacy_nonzero = legacy_nonzero,
      normalized_nonzero = normalized_nonzero,
      filter_decision_equal = legacy_nonzero == normalized_nonzero,
      original_na = original_na,
      original_sentinel = original_sentinel
    )
  readr::write_tsv(normalization_masks,
                   file.path(output_dir, "missing_normalization_masks.tsv"))
  message("Compress and index the normalized TSV")
  bgz_path <- Rsamtools::bgzip(tsv_path, dest = paste0(tsv_path, ".bgz"),
                               overwrite = FALSE)
  Rsamtools::indexTabix(bgz_path, seq = 1L, start = 2L, end = 2L, skip = 1L,
                       zeroBased = FALSE)
  unlink(tsv_path)
  copied <- file.copy(file.path(input_dir, copy_names), file.path(output_dir, copy_names),
                       overwrite = FALSE, copy.mode = TRUE)
  if (!all(copied)) stop("Could not copy all source fixture files", call. = FALSE)
  summary <- list(
    source_input_dir = input_dir,
    samples = n_samples,
    tsv_variants = n_variants,
    pgen_variants = n_variants + 4L,
    query_chromosome = "1",
    query_start = min(positions),
    query_end = max(positions),
    first_chr1_pgen_index = 3L,
    last_chr1_pgen_index = n_variants + 2L,
    missing_convention = "NA; original -1 values were normalized before filtering",
    original_na_cells = sum(original_na),
    original_sentinel_cells = sum(original_sentinel),
    normalized_hardcall_counts = as.list(hardcall_counts),
    observed_values_equal = observed_values_equal,
    missing_masks_equal = missing_masks_equal,
    pgen_roundtrip_exact = TRUE,
    pgen_roundtrip_variants_checked = roundtrip_variants_checked,
    pgen_roundtrip_cells_checked = roundtrip_cells_checked,
    legacy_variants_retained = sum(legacy_nonzero),
    normalized_variants_retained = sum(normalized_nonzero),
    filter_masks_equal = identical(legacy_nonzero, normalized_nonzero),
    filter_decisions_changed = sum(legacy_nonzero != normalized_nonzero),
    filter_expression = "rowSums(round(genotype, 0), na.rm = TRUE) != 0; evaluated by row",
    validation_scope = paste("Normalization and exact PGEN roundtrip checks cover every chr1 cell.",
                             "Separate benchmark checks compare both complete reader outputs."),
    flanks = "PGEN/PVAR only: two chr2 records before chr1 and two after chr1",
    original_fixture_json = "fixture.json is an unchanged copy describing the source fixture"
  )
  jsonlite::write_json(summary, file.path(output_dir, "pgen_preparation.json"),
                       auto_unbox = TRUE, pretty = TRUE)
  message("Prepared ", n_samples, " samples, ", n_variants,
          " matched chr1 variants; changed filter decisions: ",
          summary$filter_decisions_changed)
  invisible(summary)
}

if (sys.nframe() == 0L) {
  benchmark_library <- Sys.getenv("BENCH_LIBRARY", unset = "")
  if (nzchar(benchmark_library)) .libPaths(c(benchmark_library, .libPaths()))
  arguments <- commandArgs(trailingOnly = TRUE)
  if (length(arguments) != 2L) {
    stop("Usage: Rscript prepare_pgen_inputs.R INPUT_DIR OUTPUT_DIR", call. = FALSE)
  }
  prepare_pgen_inputs(arguments[[1L]], arguments[[2L]])
}
