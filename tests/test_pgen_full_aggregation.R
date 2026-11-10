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
    status = "CONVERGED", converged = TRUE, credible_set_count = c(1L, 0L, 1L, 0L)
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
  # Reversed file order includes an empty first file and must retain native values.
  stopifnot(identical(aggregate(rev(paths)), output))
  actual <- arrow::read_parquet(output)
  expected <- dplyr::bind_rows(tables[1:2])
  stopifnot(identical(actual[names(expected)], expected),
    all(c("batch_id", "gene_id", "outcome_id", "phenotype_id", "gene_key") %in% names(actual)),
    all(actual$molecular_trait_id == actual$outcome_id),
    identical(actual$phenotype_id, index$phenotype_id[match(actual$outcome_id, index$outcome_id)]),
    sum(is.na(actual$cs_index)) == 6L,
    identical(actual$X_column_scale_factors, expected$X_column_scale_factors),
    identical(actual$pip, expected$pip),
    identical(actual$phenotype_residual_SD, c(2, 2, 2, 0, 0, 1.5, 1.5, 1.5)))

  # Old files retain an explicit missing scale; mixing old and new does not infer one.
  legacy_tables <- lapply(tables, function(table) dplyr::select(table, -phenotype_residual_SD))
  write_tables(legacy_tables)
  aggregate()
  legacy <- arrow::read_parquet(output)
  stopifnot("phenotype_residual_SD" %in% names(legacy),
    is.double(legacy$phenotype_residual_SD), all(is.na(legacy$phenotype_residual_SD)))
  mixed_tables <- tables
  mixed_tables[[1L]] <- legacy_tables[[1L]]
  write_tables(mixed_tables)
  aggregate()
  mixed <- arrow::read_parquet(output)
  stopifnot(is.double(mixed$phenotype_residual_SD),
    identical(mixed$phenotype_residual_SD, c(rep(NA_real_, 5), 1.5, 1.5, 1.5)))
  write_tables()

  # Keep the chromosome aliases accepted by the compact aggregation path.
  for (variant_id in c("2_104_A_G", "chr2_104_A_G", "chr2_104_A_G_extra",
      "chrchr2_104_A_G", "CHRchr2_104_A_G")) {
    for (chromosome in c("2", "chr2", "chrchr2", "CHRchr2")) {
      aliased <- tables
      aliased[[1L]]$variant[[4L]] <- variant_id
      aliased[[1L]]$chromosome[[4L]] <- chromosome
      write_tables(aliased)
      aggregate()
      retained <- arrow::read_parquet(output)
      stopifnot(retained$variant[[4L]] == variant_id,
        retained$chromosome[[4L]] == chromosome)
    }
  }
  write_tables()
  aggregate()

  # Required local files include empty genes. Repeated paths and foreign names fail.
  expect_error(aggregate(paths[-3L]), "missing|exactly|gene")
  expect_error(aggregate(c(paths[1:2], paths[[2L]])), "duplicate")
  alias_dir <- file.path(root, "aliases")
  dir.create(alias_dir)
  alias <- file.path(alias_dir, basename(paths[[3L]]))
  stopifnot(file.symlink(paths[[3L]], alias))
  expect_error(aggregate(c(paths, alias)), "duplicate")
  expect_error(aggregate(c(paths[1:2], "gs://bucket/empty.parquet")), "localization")
  foreign <- file.path(root, "foreign.full_susie.parquet")
  file.copy(paths[[3L]], foreign)
  expect_error(aggregate(c(paths[1:2], foreign)), "foreign|basename|gene")
  expect_error(aggregatePgenFullSusie("gs://bucket/index.tsv", paths, output), "localization")
  expect_error(aggregate(batch_size = 0L), "batch_size")
  expect_error(aggregate(destination = "gs://bucket/full.parquet"), "local")

  check_bad_table <- function(table, pattern, gene = 1L) {
    broken <- tables
    broken[[gene]] <- table
    write_tables(broken)
    expect_error(aggregate(), pattern)
    write_tables()
  }
  # Scales must be numeric and valid even on empty tables. One outcome can span batches.
  check_bad_table(dplyr::mutate(tables[[1L]], phenotype_residual_SD = "2"), "phenotype_residual_SD|numeric|schema")
  check_bad_table(dplyr::mutate(tables[[1L]], phenotype_residual_SD = Inf), "phenotype_residual_SD|finite")
  check_bad_table(dplyr::mutate(tables[[1L]], phenotype_residual_SD = -1), "phenotype_residual_SD|negative|nonnegative")
  check_bad_table(dplyr::mutate(tables[[3L]], phenotype_residual_SD = character()), "phenotype_residual_SD|numeric|schema", 3L)
  # Row three is in the next batch; preserve the previous complete output on failure.
  previous_sd_output <- readBin(output, "raw", n = file.info(output)$size)
  for (bad_scale in c(3, NA_real_)) {
    broken <- tables
    broken[[1L]]$phenotype_residual_SD[[3L]] <- bad_scale
    write_tables(broken)
    expect_error(aggregate(batch_size = 2L), "phenotype_residual_SD|inconsistent")
    stopifnot(identical(readBin(output, "raw", n = file.info(output)$size), previous_sd_output),
      identical(list.files(dirname(output)), "full_susie.parquet"))
  }
  missing_tables <- lapply(tables, function(table) dplyr::mutate(table, phenotype_residual_SD = rep(NA_real_, nrow(table))))
  write_tables(missing_tables)
  aggregate()
  stopifnot(all(is.na(arrow::read_parquet(output)$phenotype_residual_SD)))
  write_tables()
  aggregate()
  check_bad_table(dplyr::select(tables[[1L]], -X_column_scale_factors), "schema|column")
  check_bad_table(dplyr::mutate(tables[[2L]], posterior_mean = as.character(posterior_mean)), "schema|type|numeric", 2L)
  check_bad_table(dplyr::mutate(tables[[1L]], molecular_trait_id = "outcome3"), "ownership|foreign|gene")
  check_bad_table(dplyr::mutate(tables[[1L]], chromosome = "2"), "chromosome")
  check_bad_table(dplyr::mutate(tables[[1L]], alt = "C"), "alt|allele")
  check_bad_table(dplyr::mutate(tables[[1L]], ref = "T"), "ref|allele")
  check_bad_table(dplyr::mutate(tables[[1L]], X_column_scale_factors = 0), "scale")
  check_bad_table(dplyr::mutate(tables[[1L]], X_column_scale_factors = Inf), "scale")
  check_bad_table(dplyr::bind_rows(tables[[1L]], tables[[1L]][1L, ]), "duplicate")
  check_bad_table(tables[[1L]][0L, ], "missing|credible|rows")
  check_bad_table(dplyr::mutate(tables[[3L]], pip = character()), "schema|type|numeric", 3L)

  # A late row error must not replace an earlier complete output or leave a partial file.
  previous <- readBin(output, "raw", n = file.info(output)$size)
  broken <- tables
  broken[[2L]]$X_column_scale_factors[[3L]] <- NA_real_
  write_tables(broken)
  expect_error(aggregate(), "scale")
  stopifnot(identical(readBin(output, "raw", n = file.info(output)$size), previous),
    identical(list.files(dirname(output)), "full_susie.parquet"))
  write_tables()
  bad_index <- dplyr::mutate(index, gene_key = "gene-incorrect")
  write_index(bad_index)
  expect_error(aggregate(), "gene_key|hash")
  write_index(dplyr::bind_rows(index, index[1L, ]))
  expect_error(aggregate(), "duplicate")
  write_index(dplyr::mutate(index, status = "FAILED"))
  expect_error(aggregate(), "status|FAILED")

  # Empty schema differences are accepted only on the known legacy columns.
  write_index(dplyr::mutate(index, credible_set_count = 0L))
  empties <- rep(list(dplyr::mutate(InitEmptyVariantDf(), phenotype_residual_SD = double())), length(paths))
  write_tables(empties)
  aggregate()
  empty <- arrow::read_parquet(output)
  stopifnot(nrow(empty) == 0L,
    identical(names(empty), c(names(empties[[1L]]),
      "batch_id", "gene_id", "outcome_id", "phenotype_id", "gene_key")),
    is.double(empty$phenotype_residual_SD))
  # Entirely empty legacy files still emit the optional numeric scale column.
  write_tables(lapply(empties, function(table) dplyr::select(table, -phenotype_residual_SD)))
  aggregate()
  empty_legacy <- arrow::read_parquet(output)
  stopifnot(nrow(empty_legacy) == 0L, is.double(empty_legacy$phenotype_residual_SD))
  # An integer optional column is numeric and is normalized without loosening other types.
  mixed_empty <- empties
  mixed_empty[[2L]]$phenotype_residual_SD <- integer()
  mixed_empty[[3L]] <- dplyr::select(mixed_empty[[3L]], -phenotype_residual_SD)
  write_tables(mixed_empty)
  aggregate()
  stopifnot(is.double(arrow::read_parquet(output)$phenotype_residual_SD))

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
