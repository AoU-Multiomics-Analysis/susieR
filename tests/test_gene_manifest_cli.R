# Exercise the installed aggregate CLI without fitting or opening gene artifacts.
run_gene_manifest_cli_tests <- function(repo_dir = ".") {
  suppressPackageStartupMessages(library(dplyr))
  functions_dir <- Sys.getenv("SUSIER_FUNCTIONS_PATH",
                              unset = file.path(repo_dir, "R", "utils"))
  scripts_dir <- Sys.getenv("SUSIER_SCRIPTS_PATH",
                            unset = file.path(repo_dir, "R", "scripts"))
  source(file.path(functions_dir, "InitFunctions.R"))
  previous <- Sys.getenv(c("SUSIER_FUNCTIONS_PATH", "SUSIER_SCRIPTS_PATH"), unset = NA)
  Sys.setenv(SUSIER_FUNCTIONS_PATH = functions_dir, SUSIER_SCRIPTS_PATH = scripts_dir)
  on.exit({
    for (name in names(previous)) {
      if (is.na(previous[[name]])) Sys.unsetenv(name)
      else do.call(Sys.setenv, setNames(list(previous[[name]]), name))
    }
  }, add = TRUE)
  root <- tempfile("gene manifest CLI 'quoted' $literal ")
  dir.create(root)
  passed <- FALSE
  on.exit(if (passed) unlink(root, recursive = TRUE) else
    message("Gene manifest CLI fixture retained at: ", root), add = TRUE)
  gene_ids <- c("gene A/one", "gene A/one", "gene B", "gene C")
  chromosome <- c("chr1", "chr1", "1", "chr2")
  expected <- tibble::tibble(
    batch_id = c("batch_1", "batch_1", "batch_1", "batch_2"),
    chromosome = chromosome, gene_id = gene_ids,
    outcome_id = paste0("outcome", 1:4), phenotype_id = paste0("outcome", 1:4),
    phenotype_pos = c(100, 110, 200, 300),
    requested_window_start = c(90, 100, 190, 290),
    requested_window_end = c(110, 120, 210, 310),
    window_start = c(90, 100, 190, 290), window_end = c(110, 120, 210, 310),
    region = paste0(chromosome, ":", c(90, 100, 190, 290), "-", c(110, 120, 210, 310)))
  key <- function(id) paste0("gene-", digest::digest(charToRaw(enc2utf8(id)),
    algo = "sha256", serialize = FALSE))
  index <- expected |>
    dplyr::select(batch_id, chromosome, gene_id, outcome_id, phenotype_id) |>
    dplyr::mutate(gene_key = vapply(gene_id, key, character(1)), seed = 1L,
      status = c("CONVERGED", "NONCONVERGED", "CONVERGED", "CONVERGED"),
      converged = status == "CONVERGED", credible_set_count = 0L,
      reason = dplyr::if_else(converged, "", "SuSiE did not converge"))
  expected_file <- file.path(root, "expected.tsv")
  readr::write_tsv(expected, expected_file)
  own <- expected |>
    dplyr::select(batch_id, gene_id, outcome_id, phenotype_id) |>
    dplyr::slice(0L)
  summary <- dplyr::bind_cols(own, InitEmptyInCSVariantDf(TRUE))
  sets <- dplyr::bind_cols(own, InitEmptyCS())
  diagnostics <- dplyr::bind_cols(own, tibble::tibble(
    molecular_trait_id = character(), variant = character(), input_scale = double(),
    residual_sd = double(), residual_sd_ratio = double()))
  batch_paths <- lapply(seq_len(2L), function(i) {
    prefix <- file.path(root, paste0("batch", i))
    paths <- c(summary = paste0(prefix, ".parquet"), sets = paste0(prefix, ".cs.parquet"),
      index = paste0(prefix, ".index.tsv"), diagnostics = paste0(prefix, ".diagnostics.tsv"))
    arrow::write_parquet(summary, paths[["summary"]])
    arrow::write_parquet(sets, paths[["sets"]])
    readr::write_tsv(dplyr::filter(index, batch_id == paste0("batch_", i)), paths[["index"]])
    readr::write_tsv(diagnostics, paths[["diagnostics"]])
    paths
  })
  write_list <- function(values, name) {
    path <- file.path(root, paste0(name, ".list"))
    writeLines(values, path)
    path
  }
  input_list <- function(role) write_list(vapply(batch_paths, `[[`, character(1), role), role)
  base_args <- c("--expected-outcomes", expected_file,
    "--summary-files", input_list("summary"), "--credible-set-files", input_list("sets"),
    "--index-files", input_list("index"), "--diagnostic-files", input_list("diagnostics"))
  invoke <- function(name, arguments, error = NULL) {
    log <- file.path(root, paste0(name, ".log"))
    code <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
      shQuote(c(file.path(scripts_dir, "aggregate_pgen_susie.R"), base_args, arguments)),
      stdout = log, stderr = log))
    if (is.null(error) && !identical(code, 0L)) {
      stop(name, " CLI failed: ", paste(tail(readLines(log, warn = FALSE), 20L),
        collapse = "\n"), call. = FALSE)
    }
    if (!is.null(error)) {
      stopifnot(code != 0L, grepl(error, paste(readLines(log, warn = FALSE), collapse = "\n"),
        fixed = TRUE))
    }
  }
  suffixes <- list(rds = "_susie.rds", parquet = c(".parquet", ".lbf_variable.parquet", ".full_susie.parquet"),
    tsv = c(".covariate_explained_variants.tsv", ".variant_position_summary.tsv", ".outcome_index.tsv"),
    log = ".log", provenance = ".provenance.json")
  keys <- vapply(unique(gene_ids), key, character(1))
  uris <- lapply(suffixes, function(suffix) rev(unlist(lapply(keys, function(gene_key) {
    paste0("gs://nonexistent-bucket/run 'quoted' $(touch PWNED)/", gene_key, suffix)
  }), use.names = FALSE)))
  metadata_args <- unlist(lapply(names(uris), function(role) {
    c(paste0("--gene-", role, "-uris-file"), write_list(uris[[role]], paste0("gene_", role)))
  }), use.names = FALSE)
  plain <- file.path(root, "plain")
  manifests <- file.path(root, "with_manifests")
  invoke("plain", c("--output-dir", plain))
  invoke("manifests", c(metadata_args, "--output-dir", manifests))
  for (name in c("susie.parquet", "credible_sets.parquet", "outcome_index.tsv", "covariate_explained_variants.tsv")) {
    stopifnot(identical(digest::digest(file = file.path(plain, name), algo = "sha256"),
      digest::digest(file = file.path(manifests, name), algo = "sha256")))
  }
  files <- list.files(file.path(manifests, "gene_manifests"), full.names = TRUE)
  stopifnot(setequal(basename(files), c("chr1.gene_outputs.tsv", "chr2.gene_outputs.tsv")))
  result <- dplyr::bind_rows(lapply(files, readr::read_tsv,
    col_types = readr::cols(.default = readr::col_character()), show_col_types = FALSE))
  stopifnot(nrow(result) == 3L, setequal(result$gene_id, unique(gene_ids)))
  for (i in seq_len(nrow(result))) {
    stopifnot(result$gene_key[[i]] == key(result$gene_id[[i]]),
      result$chromosome[[i]] == if (result$gene_id[[i]] == "gene C") "2" else "1")
  }
  actual_uris <- unlist(result[grepl("_uri$", names(result))], use.names = FALSE)
  stopifnot(length(actual_uris) == 27L, setequal(actual_uris, unlist(uris, use.names = FALSE)))
  invoke("partial", c(metadata_args[1:2], "--output-dir", file.path(root, "partial")),
    error = "all five")
  wrong_list <- metadata_args
  wrong_list[[2L]] <- "gs://nonexistent-bucket/unlocalized-list"
  invoke("unlocalized", c(wrong_list, "--output-dir", file.path(root, "bad")),
    error = "localization")
  stopifnot(!file.exists(file.path(root, "bad")), !file.exists("PWNED"))
  passed <- TRUE
  message("Gene manifest CLI tests passed: chromosome files, opaque URIs, and byte-identical scientific outputs")
  invisible(TRUE)
}
