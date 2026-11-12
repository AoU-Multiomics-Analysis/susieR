# Small Parquet fixtures exercise aggregation without fitting a model.
run_pgen_full_aggregation_tests <- function(repo_dir = ".") {
  functions_dir <- normalizePath(Sys.getenv("SUSIER_FUNCTIONS_PATH",
    unset = file.path(repo_dir, "R", "utils")))
  scripts_dir <- normalizePath(Sys.getenv("SUSIER_SCRIPTS_PATH",
    unset = file.path(repo_dir, "R", "scripts")))
  source(file.path(functions_dir, "InitFunctions.R"))
  source(file.path(functions_dir, "PgenBatchFunctions.R"))
  stopifnot("The full-SuSiE aggregation helper must exist" =
    exists("aggregatePgenFullSusie", mode = "function"))
  root <- tempfile("full aggregation 'quoted' $literal ")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  key <- function(gene) paste0("gene-", digest::digest(charToRaw(enc2utf8(gene)),
    algo = "sha256", serialize = FALSE))
  index <- tibble::tibble(
    batch_id = c("batch1", "batch1", "batch2", "batch2"),
    chromosome = "1", gene_id = c("gene A", "gene A", "gene B", "gene C"),
    outcome_id = paste0("outcome", 1:4), phenotype_id = paste0("phenotype", 1:4),
    status = c("NONCONVERGED", "CONVERGED", "CONVERGED", "NONCONVERGED"),
    converged = c(FALSE, TRUE, TRUE, FALSE), credible_set_count = c(1L, 0L, 1L, 0L)
  ) |>
    dplyr::mutate(gene_key = vapply(gene_id, key, character(1)))
  index_file <- file.path(root, "outcome index.tsv")
  write_index <- function(value = index) readr::write_tsv(value, index_file)
  write_index()
  paths <- file.path(root, paste0(unique(index$gene_key), ".full_susie.parquet"))
  make_table <- function(outcomes, positions) {
    n <- length(outcomes)
    table <- InitEmptyVariantDf()[rep(NA_integer_, n), ]
    for (name in names(table)[vapply(table, is.numeric, logical(1))]) table[[name]] <- rep(0.25, n)
    table |>
      dplyr::mutate(molecular_trait_id = outcomes,
        variant = paste0("chr1_", positions, "_A_G"), chromosome = "1",
        position = as.character(positions), ref = "A", alt = "G",
        cs_id = dplyr::if_else(dplyr::row_number() == 1L, paste0(outcomes, "_L1"), NA_character_),
        cs_index = dplyr::if_else(dplyr::row_number() == 1L, "L1", NA_character_),
        low_purity = dplyr::if_else(dplyr::row_number() == 1L, FALSE, NA),
        region = "chr1:1-1000", pip = seq_len(n) / (n + 1),
        X_column_scale_factors = seq_len(n) / 10,
        source_marker = paste0(outcomes, ":", positions),
        phenotype_residual_SD = unname(c(outcome1 = 2, outcome2 = 0, outcome3 = 1.5)[outcomes]))
  }
  tables <- list(make_table(c("outcome1", "outcome1", "outcome1", "outcome2", "outcome2"), 101:105),
    make_table(rep("outcome3", 3), 201:203),
    dplyr::mutate(InitEmptyVariantDf(), source_marker = character(), phenotype_residual_SD = double()))
  # A gene-owned additional feature can come from a different chromosome.
  tables[[1L]]$chromosome[[4L]] <- "2"
  tables[[1L]]$variant[[4L]] <- "chr2_104_A_G"
  write_tables <- function(values = tables) {
    for (i in seq_along(paths)) arrow::write_parquet(values[[i]], paths[[i]], chunk_size = 2L)
  }
  write_tables()
  output <- file.path(root, "aggregate", "full_susie.parquet")
  aggregate <- function(files = paths, destination = output, batch_size = 2L) {
    aggregatePgenFullSusie(index_file, files, destination, batch_size = batch_size)
  }
  expect_error <- function(expr, pattern) {
    result <- tryCatch({ force(expr); NULL }, error = identity)
    stopifnot(inherits(result, "error"), grepl(pattern, conditionMessage(result), ignore.case = TRUE))
    invisible(result)
  }
  capture_warnings <- function(expr) {
    warnings <- character()
    withCallingHandlers(force(expr), warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    })
    warnings
  }
  # Keep the supplied file order, including an empty first file.
  stopifnot(identical(aggregate(rev(paths)), output))
  actual <- arrow::read_parquet(output)
  expected <- dplyr::bind_rows(tables[c(2L, 1L)])
  forward_expected <- dplyr::bind_rows(tables[1:2])
  stopifnot(identical(actual[names(expected)], expected),
    all(c("batch_id", "gene_id", "outcome_id", "phenotype_id", "gene_key") %in% names(actual)),
    all(actual$molecular_trait_id == actual$outcome_id),
    identical(actual$phenotype_id, index$phenotype_id[match(actual$outcome_id, index$outcome_id)]),
    sum(is.na(actual$cs_index)) == 6L,
    identical(actual$phenotype_residual_SD, c(1.5, 1.5, 1.5, 2, 2, 2, 0, 0)),
    identical(actual$feature_converged, c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE)))
  aggregate()
  forward <- arrow::read_parquet(output)
  stopifnot(identical(forward[names(forward_expected)], forward_expected))

  # Equal labels and differing values are ordinary rows, within or across batches.
  differing_row <- tables[[1L]][1L, ] |>
    dplyr::mutate(pip = .9, z = -3, posterior_mean = -.4, posterior_sd = .6,
                  alpha1 = .8, mu_1 = -.7, mu2_1 = .9, source_marker = "different")
  repeated_cases <- list(
    identical_across_batches = dplyr::bind_rows(tables[[1L]], tables[[1L]][1L, ]),
    differing_across_batches = dplyr::bind_rows(tables[[1L]], differing_row),
    repeated_within_batch = dplyr::bind_rows(tables[[1L]][1L, ], differing_row,
                                            tables[[1L]][-1L, ])
  )
  for (case in names(repeated_cases)) {
    repeated_tables <- tables
    repeated_tables[[1L]] <- repeated_cases[[case]]
    write_tables(repeated_tables)
    stopifnot(length(capture_warnings(aggregate(batch_size = 2L))) == 0L)
    retained <- arrow::read_parquet(output)
    expected_repeated <- dplyr::bind_rows(repeated_tables[1:2])
    stopifnot(identical(retained[names(expected_repeated)], expected_repeated))
  }
  write_tables()

  # Repeated paths, subsets, aliases, and arbitrary names retain their rows.
  aggregate(c(paths[[2L]], paths[[1L]], paths[[2L]], paths[[3L]]))
  repeated_paths <- dplyr::bind_rows(tables[c(2L, 1L, 2L)])
  stopifnot(identical(arrow::read_parquet(output)[names(repeated_paths)], repeated_paths))
  alias_dir <- file.path(root, "aliases")
  dir.create(alias_dir)
  alias <- file.path(alias_dir, basename(paths[[1L]]))
  stopifnot(file.symlink(paths[[1L]], alias))
  aggregate(c(paths[[1L]], alias))
  duplicated_source <- dplyr::bind_rows(tables[c(1L, 1L)])
  stopifnot(identical(arrow::read_parquet(output)[names(duplicated_source)], duplicated_source))
  foreign <- file.path(root, "arbitrary subset.parquet")
  stopifnot(file.copy(paths[[1L]], foreign))
  aggregate(foreign)
  stopifnot(identical(arrow::read_parquet(output)[names(tables[[1L]])], tables[[1L]]))
  aggregate(paths[[3L]])
  stopifnot(nrow(arrow::read_parquet(output)) == 0L)

  # Supplied convergence values stay intact; only missing values use the index.
  convergence_tables <- tables
  convergence_tables[[1L]]$feature_converged <- c(NA, TRUE, FALSE, FALSE, NA)
  convergence_tables[[3L]]$feature_converged <- logical()
  write_tables(convergence_tables)
  aggregate()
  mixed_convergence <- arrow::read_parquet(output)
  stopifnot(identical(mixed_convergence$feature_converged,
    c(FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)),
    identical(mixed_convergence[names(forward_expected)], forward_expected))

  # Missing legacy scale columns become numeric NA. Existing values are copied.
  legacy_tables <- lapply(tables, function(table) dplyr::select(table, -phenotype_residual_SD))
  write_tables(legacy_tables)
  aggregate()
  legacy <- arrow::read_parquet(output)
  stopifnot(is.double(legacy$phenotype_residual_SD), all(is.na(legacy$phenotype_residual_SD)))
  mixed_tables <- tables
  mixed_tables[[1L]] <- legacy_tables[[1L]]
  write_tables(mixed_tables)
  aggregate()
  mixed <- arrow::read_parquet(output)
  stopifnot(is.double(mixed$phenotype_residual_SD),
    identical(mixed$phenotype_residual_SD, c(rep(NA_real_, 5), 1.5, 1.5, 1.5)))

  # An absent optional column must use the supplied column type in other files.
  character_scale <- tables
  character_scale[[1L]] <- dplyr::select(character_scale[[1L]], -phenotype_residual_SD)
  character_scale[[2L]]$phenotype_residual_SD <- as.character(character_scale[[2L]]$phenotype_residual_SD)
  character_scale[[3L]] <- dplyr::select(character_scale[[3L]], -phenotype_residual_SD)
  write_tables(character_scale)
  aggregate()
  character_merged <- arrow::read_parquet(output)
  character_expected <- dplyr::bind_rows(character_scale[1:2])
  stopifnot(identical(character_merged[names(character_expected)], character_expected),
    identical(character_merged$phenotype_residual_SD, c(rep(NA_character_, 5), rep("1.5", 3L))))

  # Use the column union. Missing native or extra columns receive NA values.
  union_tables <- tables
  union_tables[[1L]] <- dplyr::select(union_tables[[1L]], -alpha2, -X_column_scale_factors)
  union_tables[[2L]] <- union_tables[[2L]] |>
    dplyr::select(-source_marker) |>
    dplyr::mutate(extra_annotation = c("first", "second", "third"))
  write_tables(union_tables)
  aggregate()
  union_expected <- dplyr::bind_rows(union_tables[1:2])
  stopifnot(identical(arrow::read_parquet(output)[names(union_expected)], union_expected))
  arbitrary <- tibble::tibble(source_note = c("first", "second"), pip = c(-2, 3))
  arbitrary_path <- file.path(root, "table without native fields.parquet")
  arrow::write_parquet(arbitrary, arbitrary_path)
  aggregate(arbitrary_path)
  arbitrary_merged <- arrow::read_parquet(output)
  stopifnot(identical(arbitrary_merged[names(arbitrary)], arbitrary),
    all(is.na(arbitrary_merged$gene_id)), all(is.na(arbitrary_merged$feature_converged)))

  # Scientific inconsistencies do not prevent a merge or change source values.
  inconsistent <- tables
  inconsistent[[1L]] <- inconsistent[[1L]] |>
    dplyr::mutate(molecular_trait_id = c("outcome1", "outcome3", "unknown", NA, ""),
      variant = c("not a variant", "chrX_999_T_C", "same", "same", "same"),
      chromosome = "different", position = "not a coordinate", ref = "T", alt = "C",
      pip = c(-1, 2, NA_real_, Inf, .25), X_column_scale_factors = c(0, -1, Inf, NA_real_, 1),
      phenotype_residual_SD = c(2, -1, Inf, NA_real_, 3),
      cs_id = "unrelated CS", cs_index = "not native", low_purity = FALSE)
  permissive_index <- index |>
    dplyr::mutate(status = "FAILED", credible_set_count = -1L,
      chromosome = "X", gene_key = "supplied key")
  write_index(permissive_index)
  write_tables(inconsistent)
  stopifnot(length(capture_warnings(aggregate())) == 0L)
  retained <- arrow::read_parquet(output)
  inconsistent_expected <- dplyr::bind_rows(inconsistent[1:2])
  stopifnot(identical(retained[names(inconsistent_expected)], inconsistent_expected),
    identical(retained$gene_id[1:2], c("gene A", "gene B")),
    all(is.na(retained$gene_id[3:5])), all(is.na(retained$feature_converged[3:5])))

  # Conflicting index duplicates cannot multiply rows or supply guessed metadata.
  ambiguous_index <- dplyr::bind_rows(index, index[1L, ] |>
    dplyr::mutate(gene_id = "other gene", phenotype_id = "other phenotype", converged = TRUE))
  write_index(ambiguous_index)
  write_tables()
  stopifnot(length(capture_warnings(aggregate())) == 0L)
  ambiguous <- arrow::read_parquet(output)
  stopifnot(identical(ambiguous[names(forward_expected)], forward_expected),
    all(is.na(ambiguous$feature_converged[1:3])), all(is.na(ambiguous$gene_id[1:3])),
    identical(ambiguous$feature_converged[4:8], rep(TRUE, 5L)))

  # Preserve source ownership and metadata, including supplied NA values.
  supplied <- tables
  supplied[[1L]] <- supplied[[1L]] |>
    dplyr::mutate(batch_id = c("source batch", NA, "source batch", "source batch", "source batch"),
      gene_id = "source gene", outcome_id = "source outcome", phenotype_id = "source phenotype",
      gene_key = "source key", status = "source status", converged = FALSE,
      credible_set_count = 99L, extra_metadata = "keep this")
  write_index()
  write_tables(supplied)
  aggregate()
  retained <- arrow::read_parquet(output)
  stopifnot(identical(retained[1:5, names(supplied[[1L]])], supplied[[1L]]),
    identical(retained$batch_id[6:8], rep("batch2", 3L)),
    all(is.na(retained$extra_metadata[6:8])))

  # A sparse index can provide convergence without ownership metadata.
  write_index(dplyr::select(index, outcome_id, converged))
  write_tables()
  aggregate()
  sparse <- arrow::read_parquet(output)
  stopifnot(identical(sparse[names(forward_expected)], forward_expected),
    all(is.na(sparse$gene_id)), identical(sparse$outcome_id, sparse$molecular_trait_id),
    identical(sparse$feature_converged, c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE)))
  write_index(index[0L, ])
  aggregate()
  no_index_rows <- arrow::read_parquet(output)
  stopifnot(identical(no_index_rows[names(forward_expected)], forward_expected),
    all(is.na(no_index_rows$gene_id)), all(is.na(no_index_rows$feature_converged)))
  write_index()
  write_tables()
  aggregate()

  # I/O failures preserve the complete output and remove any temporary output.
  previous <- readBin(output, "raw", n = file.info(output)$size)
  check_unchanged_output <- function() {
    stopifnot(identical(readBin(output, "raw", n = file.info(output)$size), previous),
      identical(list.files(dirname(output), all.files = TRUE, no.. = TRUE), "full_susie.parquet"))
  }
  expect_error(aggregate(c(paths[[1L]], file.path(root, "missing later.parquet"))), "local|missing|read|exist")
  check_unchanged_output()
  corrupt <- file.path(root, "corrupt later.parquet")
  bytes <- readBin(paths[[2L]], "raw", n = file.info(paths[[2L]])$size)
  bytes[5:20] <- as.raw(255L)
  writeBin(bytes, corrupt)
  expect_error(aggregate(c(paths[[1L]], corrupt)), ".")
  check_unchanged_output()
  expect_error(aggregate(c(paths[[1L]], "gs://bucket/later.parquet")), "localization")
  expect_error(aggregatePgenFullSusie("gs://bucket/index.tsv", paths, output), "localization")
  for (size in list(0L, 1.5, Inf, NA_real_, "2")) expect_error(aggregate(batch_size = size), "batch_size")
  expect_error(aggregate(destination = "gs://bucket/full.parquet"), "local")
  expect_error(aggregate(destination = paths[[1L]]), "replace|overwrite|input")
  expect_error(aggregate(destination = index_file), "replace|overwrite|input")
  check_unchanged_output()

  # Empty native schemas remain compatible with old and new optional columns.
  empties <- rep(list(dplyr::mutate(InitEmptyVariantDf(), phenotype_residual_SD = double())), length(paths))
  write_tables(empties)
  aggregate()
  empty <- arrow::read_parquet(output)
  stopifnot(nrow(empty) == 0L, all(names(empties[[1L]]) %in% names(empty)),
    is.double(empty$phenotype_residual_SD), is.logical(empty$feature_converged))
  mixed_empty <- empties
  mixed_empty[[1L]] <- dplyr::select(mixed_empty[[1L]], -phenotype_residual_SD)
  mixed_empty[[2L]]$phenotype_residual_SD <- integer()
  mixed_empty[[2L]]$position <- character()
  mixed_empty[[2L]]$low_purity <- logical()
  mixed_empty[[3L]]$extra_empty_metadata <- character()
  write_tables(mixed_empty)
  aggregate()
  empty_union <- arrow::read_parquet(output)
  stopifnot(nrow(empty_union) == 0L, "extra_empty_metadata" %in% names(empty_union),
    is.numeric(empty_union$phenotype_residual_SD), is.logical(empty_union$feature_converged))
  write_tables(empties)

  # The CLI must treat cloud output locations as metadata while opening local Files.
  write_index(dplyr::mutate(index, credible_set_count = 0L))
  expected_file <- file.path(root, "expected outcomes.tsv")
  expected_plan <- index |>
    dplyr::mutate(phenotype_pos = 100, requested_window_start = 1,
      requested_window_end = 1000, window_start = 1, window_end = 1000,
      region = "chr1:1-1000")
  readr::write_tsv(expected_plan, expected_file)
  owners <- dplyr::select(index, batch_id, gene_id, outcome_id, phenotype_id)[0L, ]
  compact_paths <- lapply(unique(index$batch_id), function(batch) {
    prefix <- file.path(root, batch)
    files <- c(summary = paste0(prefix, ".parquet"), sets = paste0(prefix, ".cs.parquet"),
      index = paste0(prefix, ".index.tsv"), diagnostics = paste0(prefix, ".diagnostics.tsv"))
    arrow::write_parquet(dplyr::bind_cols(owners, InitEmptyInCSVariantDf(TRUE)), files[["summary"]])
    arrow::write_parquet(dplyr::bind_cols(owners, InitEmptyCS()), files[["sets"]])
    readr::write_tsv(dplyr::filter(index, batch_id == batch) |>
      dplyr::mutate(credible_set_count = 0L, seed = 1L, reason = ""), files[["index"]])
    readr::write_tsv(dplyr::bind_cols(owners, tibble::tibble(molecular_trait_id = character(),
      variant = character(), input_scale = double(), residual_sd = double(),
      residual_sd_ratio = double())), files[["diagnostics"]])
    files
  })
  list_file <- function(values, name) {
    path <- file.path(root, paste0(name, " list.txt"))
    writeLines(values, path)
    path
  }
  args <- c("--expected-outcomes", expected_file,
    "--summary-files", list_file(vapply(compact_paths, `[[`, character(1), "summary"), "summary"),
    "--credible-set-files", list_file(vapply(compact_paths, `[[`, character(1), "sets"), "sets"),
    "--index-files", list_file(vapply(compact_paths, `[[`, character(1), "index"), "index"),
    "--diagnostic-files", list_file(vapply(compact_paths, `[[`, character(1), "diagnostics"), "diagnostics"),
    "--full-susie-files", list_file(paths, "full"),
    "--output-dir", file.path(root, "cli aggregate"))
  suffixes <- list(rds = "_susie.rds",
    parquet = c(".parquet", ".lbf_variable.parquet", ".full_susie.parquet"),
    tsv = c(".covariate_explained_variants.tsv", ".variant_position_summary.tsv", ".outcome_index.tsv"),
    log = ".log", provenance = ".provenance.json")
  cloud_prefix <- "gs://nonexistent-bucket/run 'quoted' $literal/"
  metadata_args <- unlist(lapply(names(suffixes), function(role) {
    uris <- unlist(lapply(suffixes[[role]], function(suffix)
      paste0(cloud_prefix, unique(index$gene_key), suffix)), use.names = FALSE)
    c(paste0("--gene-", role, "-uris-file"), list_file(uris, paste0("metadata ", role)))
  }), use.names = FALSE)
  previous_functions <- Sys.getenv("SUSIER_FUNCTIONS_PATH", unset = NA_character_)
  Sys.setenv(SUSIER_FUNCTIONS_PATH = functions_dir)
  on.exit(if (is.na(previous_functions)) Sys.unsetenv("SUSIER_FUNCTIONS_PATH") else
    Sys.setenv(SUSIER_FUNCTIONS_PATH = previous_functions), add = TRUE)
  log <- file.path(root, "cli.log")
  code <- system2(file.path(R.home("bin"), "Rscript"),
    shQuote(c(file.path(scripts_dir, "aggregate_pgen_susie.R"), args, metadata_args)), stdout = log, stderr = log)
  if (!identical(code, 0L)) stop(paste(readLines(log, warn = FALSE), collapse = "\n"))
  stopifnot(nrow(arrow::read_parquet(file.path(root, "cli aggregate", "full_susie.parquet"))) == 0L)
  manifest <- readr::read_tsv(file.path(root, "cli aggregate", "gene_manifests", "chr1.gene_outputs.tsv"),
    show_col_types = FALSE)
  stopifnot(identical(manifest$full_parquet_uri,
    paste0(cloud_prefix, unique(index$gene_key), ".full_susie.parquet")))

  # The normal CLI merges supplied rows despite failed/duplicate index entries
  # and unrelated expected-input content. Output URI lists remain metadata.
  writeLines("readable placeholder", expected_file)
  cli_summaries <- lapply(seq_along(compact_paths), function(i) {
    batch <- unique(index$batch_id)[[i]]
    own <- dplyr::filter(index, batch_id == batch)
    failed <- dplyr::mutate(own, status = "FAILED", converged = NA,
      credible_set_count = -1L, seed = 1L, reason = "supplied status")
    readr::write_tsv(dplyr::bind_rows(failed, failed[1L, ]), compact_paths[[i]][["index"]])
    row <- InitEmptyInCSVariantDf(TRUE)[NA_integer_, ] |>
      dplyr::mutate(molecular_trait_id = own$outcome_id[[1L]],
        variant = "supplied label", pip = -2, source_marker = paste0("source", i))
    summary <- dplyr::bind_cols(
      dplyr::select(own[1L, ], batch_id, gene_id, outcome_id, phenotype_id), row)
    arrow::write_parquet(summary, compact_paths[[i]][["summary"]])
    summary
  })
  pgen_values <- list("pgen-chromosomes-file" = "2",
    "pgen-uris-file" = paste0(cloud_prefix, "chr2.pgen"),
    "pvar-uris-file" = paste0(cloud_prefix, "chr2.pvar"),
    "psam-uris-file" = paste0(cloud_prefix, "chr2.psam"))
  pgen_args <- unlist(lapply(names(pgen_values), function(name) {
    c(paste0("--", name), list_file(pgen_values[[name]], name))
  }), use.names = FALSE)
  relaxed_args <- args
  relaxed_dir <- file.path(root, "cli supplied metadata")
  relaxed_args[[match("--output-dir", relaxed_args) + 1L]] <- relaxed_dir
  code <- system2(file.path(R.home("bin"), "Rscript"),
    shQuote(c(file.path(scripts_dir, "aggregate_pgen_susie.R"),
      relaxed_args, metadata_args, pgen_args)), stdout = log, stderr = log)
  if (!identical(code, 0L)) stop(paste(readLines(log, warn = FALSE), collapse = "\n"))
  cli_expected <- dplyr::bind_rows(cli_summaries)
  cli_actual <- arrow::read_parquet(file.path(relaxed_dir, "susie.parquet"))
  cli_index <- readr::read_tsv(file.path(relaxed_dir, "outcome_index.tsv"),
    show_col_types = FALSE)
  pgen_manifest <- readr::read_tsv(file.path(relaxed_dir, "pgen_manifest.tsv"),
    col_types = "cccc", show_col_types = FALSE)
  stopifnot(identical(cli_actual[names(cli_expected)], cli_expected),
    nrow(cli_index) == 6L, all(cli_index$status == "FAILED"),
    file.exists(file.path(relaxed_dir, "gene_manifests", "chr1.gene_outputs.tsv")),
    identical(pgen_manifest$chromosome, "2"),
    identical(pgen_manifest$pgen, pgen_values[["pgen-uris-file"]]))

  # Full-only mode consumes existing files directly, without compact input tables.
  write_index()
  write_tables()
  full_only <- c("--full-only", "--outcome-index", index_file,
    "--full-susie-files", list_file(rev(paths), "existing full"),
    "--output-dir", file.path(root, "full only"))
  code <- system2(file.path(R.home("bin"), "Rscript"),
    shQuote(c(file.path(scripts_dir, "aggregate_pgen_susie.R"), full_only)), stdout = log, stderr = log)
  if (!identical(code, 0L)) stop(paste(readLines(log, warn = FALSE), collapse = "\n"))
  stopifnot(identical(arrow::read_parquet(file.path(root, "full only", "full_susie.parquet"))[
    names(expected)], expected))
  code <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    shQuote(c(file.path(scripts_dir, "aggregate_pgen_susie.R"), full_only,
      "--expected-outcomes", expected_file)), stdout = log, stderr = log))
  stopifnot(code != 0L, any(grepl("cannot combine", readLines(log, warn = FALSE), fixed = TRUE)))
  code <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    shQuote(c(file.path(scripts_dir, "aggregate_pgen_susie.R"),
      "--outcome-index", index_file, "--output-dir", file.path(root, "invalid mode"))),
    stdout = log, stderr = log))
  stopifnot(code != 0L, any(grepl("requires --full-only", readLines(log, warn = FALSE), fixed = TRUE)))
  code <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    shQuote(c(file.path(scripts_dir, "aggregate_pgen_susie.R"), full_only,
      metadata_args[1:2])), stdout = log, stderr = log))
  stopifnot(code != 0L, any(grepl("cannot combine", readLines(log, warn = FALSE), fixed = TRUE)))
  message("Full-SuSiE aggregation tests passed")
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  run_pgen_full_aggregation_tests(if (length(args)) args[[1L]] else ".")
}
