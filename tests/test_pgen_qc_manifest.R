# Shared QC must group complete triplets, and reusable manifests must keep URIs as data.
run_pgen_qc_manifest_tests <- function(repo_dir = ".") {
  suppressPackageStartupMessages(library(dplyr))
  functions_dir <- Sys.getenv("SUSIER_FUNCTIONS_PATH",
                              unset = file.path(repo_dir, "R", "utils"))
  scripts_dir <- Sys.getenv("SUSIER_SCRIPTS_PATH",
                            unset = file.path(repo_dir, "R", "scripts"))
  source(file.path(functions_dir, "PgenBatchFunctions.R"))
  stopifnot("The planner must group unique selected PGEN triplets" =
              exists("writePgenQcPlan", mode = "function"),
            "Aggregation must write a reusable PGEN manifest" =
              exists("writePgenOutputManifest", mode = "function"))
  root <- tempfile("pgen QC metadata 'quoted' $literal ")
  dir.create(root)
  passed <- FALSE
  on.exit(if (passed) unlink(root, recursive = TRUE) else
    message("PGEN QC metadata fixture retained at: ", root), add = TRUE)
  error <- function(code, pattern) {
    result <- tryCatch({ force(code); NULL }, error = identity)
    stopifnot(inherits(result, "error"),
              grepl(pattern, conditionMessage(result), ignore.case = TRUE))
    invisible(result)
  }
  uri <- "gs://unopened-bucket/shared 'quoted' \"double\" $(touch PWNED) `touch PWNED`"
  manifest <- tibble::tibble(
    chromosome = c("10", "chr2", "5", "chr3", "6", "4"),
    pgen = paste0(c("gs://unopened-bucket/unused10", uri, "gs://unopened-bucket/b", uri,
                   "gs://unopened-bucket/unused6", uri), ".pgen"),
    pvar = paste0(c("gs://unopened-bucket/unused10", uri, "gs://unopened-bucket/b", uri,
                   "gs://unopened-bucket/unused6", paste0(uri, "-other")), ".pvar"),
    psam = rep(paste0(uri, ".psam"), 6L))
  plan_dir <- file.path(root, "unit plan")
  planned <- writePgenQcPlan(manifest, c("3", "chr5", "2", "chr3", "4"), plan_dir)
  expected_unique <- manifest[c(2L, 3L, 6L), ]
  stopifnot(identical(planned$indices, c(0L, 1L, 0L, 0L, 2L)),
    identical(lapply(planned$unique_pgens, identity), lapply(expected_unique, identity)),
    identical(readLines(planned$unique_manifest), unname(apply(as.data.frame(expected_unique), 1L, paste, collapse = "\t"))),
    identical(jsonlite::read_json(planned$batch_pgen_indices, simplifyVector = TRUE),
              c(0L, 1L, 0L, 0L, 2L)))
  stopifnot("Planner must return grouped batch indices" =
              is.character(planned$pgen_batch_indices) && is.character(planned$batch_within_pgen_indices))
  stopifnot(identical(jsonlite::read_json(planned$pgen_batch_indices, simplifyVector = TRUE),
                      list(c(0L, 2L, 3L), 1L, 4L)),
    identical(jsonlite::read_json(planned$batch_within_pgen_indices, simplifyVector = TRUE),
              c(0L, 0L, 1L, 2L, 0L)))
  error(writePgenQcPlan(manifest, c("2", "7"), file.path(root, "absent")), "chromosome|manifest")
  error(writePgenQcPlan(manifest, character(), file.path(root, "empty")), "nonempty|batch")
  different_psam <- manifest
  different_psam$pvar[[6L]] <- manifest$pvar[[2L]]
  different_psam$psam[[6L]] <- "gs://unopened-bucket/other.psam"
  psam_plan <- writePgenQcPlan(different_psam, c("2", "3", "4"), file.path(root, "different psam"))
  stopifnot(identical(psam_plan$indices, c(0L, 0L, 1L)), nrow(psam_plan$unique_pgens) == 2L)

  # Repeated chromosome aliases may refer to one exact output triplet only.
  output_file <- file.path(root, "output manifest.tsv")
  chromosomes <- c("chr2", "2", "chr5", "3", "4")
  pgens <- paste0(c(uri, uri, "gs://unopened-bucket/b", uri, uri), ".pgen")
  pvars <- paste0(c(uri, uri, "gs://unopened-bucket/b", uri, paste0(uri, "-other")), ".pvar")
  psams <- rep(paste0(uri, ".psam"), 5L)
  written <- writePgenOutputManifest(chromosomes, pgens, pvars, psams, output_file)
  result <- readPgenManifest(written)
  stopifnot(identical(result$chromosome, c("2", "5", "3", "4")),
    identical(result$pgen, pgens[c(1L, 3L, 4L, 5L)]),
    identical(result$pvar, pvars[c(1L, 3L, 4L, 5L)]),
    identical(result$psam, psams[c(1L, 3L, 4L, 5L)]),
    identical(readLines(written)[[1L]], "chromosome\tpgen\tpvar\tpsam"))
  invalid_outputs <- list(
    unequal = list(pgen_uris = pgens[-1L]),
    conflicting_pgen = list(pgen_uris = replace(pgens, 2L, "gs://bucket/conflict.pgen")),
    conflicting_pvar = list(pvar_uris = replace(pvars, 2L, "gs://bucket/conflict.pvar")),
    conflicting_psam = list(psam_uris = replace(psams, 2L, "gs://bucket/conflict.psam")),
    relative = list(pgen_uris = replace(pgens, 1:2, "relative.pgen")),
    empty_uri = list(pgen_uris = replace(pgens, 1:2, "gs://")),
    multiline = list(psam_uris = replace(psams, 1:2, "gs://bucket/a\nb.psam")),
    sex = list(chromosomes = replace(chromosomes, 1:2, "X")))
  for (name in names(invalid_outputs)) {
    arguments <- utils::modifyList(list(chromosomes = chromosomes, pgen_uris = pgens,
      pvar_uris = pvars, psam_uris = psams,
      output_file = file.path(root, paste0(name, ".tsv"))), invalid_outputs[[name]])
    error(do.call(writePgenOutputManifest, arguments), "length|triplet|location|absolute|chromosome|autosome|nonempty|single")
    stopifnot(!file.exists(arguments$output_file))
  }

  previous <- Sys.getenv(c("SUSIER_FUNCTIONS_PATH", "SUSIER_SCRIPTS_PATH"), unset = NA)
  Sys.setenv(SUSIER_FUNCTIONS_PATH = functions_dir, SUSIER_SCRIPTS_PATH = scripts_dir)
  on.exit({
    for (name in names(previous)) {
      if (is.na(previous[[name]])) Sys.unsetenv(name)
      else do.call(Sys.setenv, setNames(list(previous[[name]]), name))
    }
  }, add = TRUE)
  invoke <- function(script, name, arguments, expected_error = NULL) {
    log <- file.path(root, paste0(name, ".log"))
    status <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
      shQuote(c(file.path(scripts_dir, script), arguments)), stdout = log, stderr = log))
    detail <- paste(readLines(log, warn = FALSE), collapse = "\n")
    if (is.null(expected_error) && !identical(status, 0L)) {
      stop(name, " CLI failed: ", detail, call. = FALSE)
    }
    if (!is.null(expected_error)) {
      stopifnot(status != 0L, grepl(expected_error, detail, ignore.case = TRUE))
    }
    invisible(status)
  }
  bed <- file.path(root, "phenotypes.bed")
  qtl <- file.path(root, "qtl.tsv")
  manifest_file <- file.path(root, "input manifest.tsv")
  writeLines(c("#chr\tstart\tend\tphenotype_id\tS1",
    "chr2\t99\t101\ta\t1", "chr2\t199\t201\tb\t2", "3\t299\t301\tc\t3",
    "4\t399\t401\td\t4", "5\t499\t501\te\t5"), bed)
  readr::write_tsv(tibble::tibble(phenotype_id = letters[1:5], qval = 0.01), qtl)
  readr::write_tsv(manifest, manifest_file, quote = "all")
  cli_plan <- file.path(root, "CLI plan")
  plan_args <- c("--phenotype-bed", bed, "--qtl-results", qtl,
    "--pgen-manifest", manifest_file, "--cis-distance", "10", "--batch-size", "1")
  invoke("plan_pgen_susie.R", "plan", c(plan_args, "--output-dir", cli_plan))
  stopifnot(identical(readLines(file.path(cli_plan, "chromosome_pgens.tsv")),
                     unname(apply(as.data.frame(manifest), 1L, paste, collapse = "\t"))),
    identical(readLines(file.path(cli_plan, "unique_pgens.tsv")),
              unname(apply(as.data.frame(expected_unique), 1L, paste, collapse = "\t"))),
    identical(jsonlite::read_json(file.path(cli_plan, "batch_pgen_indices.json"), simplifyVector = TRUE),
              c(0L, 0L, 0L, 2L, 1L)),
    identical(jsonlite::read_json(file.path(cli_plan, "batch_chromosome_indices.json"), simplifyVector = TRUE),
              c(1L, 1L, 3L, 5L, 2L)))
  stopifnot(identical(jsonlite::read_json(file.path(cli_plan, "pgen_batch_indices.json"), simplifyVector = TRUE),
                      list(0:2, 4L, 3L)),
    identical(jsonlite::read_json(file.path(cli_plan, "batch_within_pgen_indices.json"), simplifyVector = TRUE),
              c(0L, 1L, 2L, 0L, 0L)))
  for (arguments in list(c("--mac=-1"), c("--mac=0.5"), c("--min-call-rate=-0.1"),
                         c("--min-call-rate=1.01"), c("--min-call-rate=NaN"))) {
    invoke("plan_pgen_susie.R", "bad threshold", c(plan_args, arguments,
      "--output-dir", file.path(root, "bad threshold")), "MAC|call.rate")
  }
  off_plan <- file.path(root, "off plan")
  invoke("plan_pgen_susie.R", "off plan", c(plan_args, "--mac=0", "--min-call-rate=0", "--output-dir", off_plan))
  stopifnot(identical(readLines(file.path(cli_plan, "unique_pgens.tsv")),
                     readLines(file.path(off_plan, "unique_pgens.tsv"))))
  chromosome_file <- file.path(root, "chromosomes.txt")
  writeLines(manifest$chromosome, chromosome_file)
  legacy_plan <- file.path(root, "legacy plan")
  invoke("plan_pgen_susie.R", "legacy plan", c("--phenotype-bed", bed, "--qtl-results", qtl,
    "--chromosomes-file", chromosome_file, "--cis-distance", "10", "--batch-size", "1",
    "--output-dir", legacy_plan))
  stopifnot(identical(readLines(file.path(cli_plan, "expected_outcomes.tsv")),
    readLines(file.path(legacy_plan, "expected_outcomes.tsv"))),
    !file.exists(file.path(legacy_plan, "unique_pgens.tsv")))

  # Run the real aggregate CLI. Output locations do not exist and must not be opened.
  source(file.path(functions_dir, "InitFunctions.R"))
  expected <- pgenReadTextTable(file.path(cli_plan, "expected_outcomes.tsv"))
  own <- expected |>
    dplyr::select(batch_id, gene_id, outcome_id, phenotype_id) |>
    dplyr::slice(0L)
  summary <- dplyr::bind_cols(own, InitEmptyInCSVariantDf(TRUE))
  sets <- dplyr::bind_cols(own, InitEmptyCS())
  diagnostics <- dplyr::bind_cols(own, tibble::tibble(molecular_trait_id = character(),
    variant = character(), input_scale = double(), residual_sd = double(), residual_sd_ratio = double()))
  indices <- expected |>
    dplyr::select(batch_id, chromosome, gene_id, outcome_id, phenotype_id) |>
    dplyr::mutate(gene_key = paste0("gene-", gene_id), seed = 1L, status = "CONVERGED",
      converged = TRUE, credible_set_count = 0L, reason = "")
  batches <- lapply(seq_len(nrow(indices)), function(i) {
    prefix <- file.path(root, paste0("batch", i))
    paths <- c(summary = paste0(prefix, ".parquet"), sets = paste0(prefix, ".cs.parquet"),
      index = paste0(prefix, ".index.tsv"), diagnostics = paste0(prefix, ".diagnostics.tsv"))
    arrow::write_parquet(summary, paths[["summary"]])
    arrow::write_parquet(sets, paths[["sets"]])
    readr::write_tsv(indices[i, ], paths[["index"]])
    readr::write_tsv(diagnostics, paths[["diagnostics"]])
    paths
  })
  write_list <- function(values, name) {
    path <- file.path(root, paste0(name, ".list"))
    writeLines(values, path)
    path
  }
  input_list <- function(role) write_list(vapply(batches, `[[`, character(1), role), role)
  aggregate_args <- c("--expected-outcomes", file.path(cli_plan, "expected_outcomes.tsv"),
    "--summary-files", input_list("summary"), "--credible-set-files", input_list("sets"),
    "--index-files", input_list("index"), "--diagnostic-files", input_list("diagnostics"))
  output_triplets <- expected_unique[c(1L, 1L, 1L, 3L, 2L), ]
  output_metadata <- list(pgen_chromosomes_file = c("chr2", "2", "chr3", "4", "5"),
    pgen_uris_file = output_triplets$pgen, pvar_uris_file = output_triplets$pvar,
    psam_uris_file = output_triplets$psam)
  metadata_paths <- mapply(write_list, output_metadata, names(output_metadata), SIMPLIFY = FALSE)
  metadata_args <- unlist(lapply(names(metadata_paths), function(name) {
    c(paste0("--", gsub("_", "-", name)), metadata_paths[[name]])
  }), use.names = FALSE)
  plain <- file.path(root, "plain aggregate")
  prepared <- file.path(root, "prepared aggregate")
  invoke("aggregate_pgen_susie.R", "plain aggregate", c(aggregate_args, "--output-dir", plain))
  invoke("aggregate_pgen_susie.R", "prepared aggregate", c(aggregate_args, metadata_args, "--output-dir", prepared))
  stopifnot(!file.exists(file.path(plain, "pgen_manifest.tsv")))
  result <- readPgenManifest(file.path(prepared, "pgen_manifest.tsv"))
  stopifnot(identical(result$chromosome, c("2", "3", "4", "5")),
    identical(result$pgen, output_triplets$pgen[c(1L, 3L, 4L, 5L)]),
    identical(result$pvar, output_triplets$pvar[c(1L, 3L, 4L, 5L)]))
  for (name in c("susie.parquet", "credible_sets.parquet", "outcome_index.tsv", "covariate_explained_variants.tsv")) {
    stopifnot(identical(digest::digest(file = file.path(plain, name), algo = "sha256"),
                       digest::digest(file = file.path(prepared, name), algo = "sha256")))
  }
  invoke("aggregate_pgen_susie.R", "partial metadata", c(aggregate_args, metadata_args[1:2],
    "--output-dir", file.path(root, "partial metadata")), "all four")
  bad_args <- metadata_args
  bad_args[[2L]] <- "gs://bucket/unlocalized.list"
  invoke("aggregate_pgen_susie.R", "unlocalized metadata", c(aggregate_args, bad_args,
    "--output-dir", file.path(root, "unlocalized metadata")), "localization")
  bad_args[[2L]] <- file.path(root, "missing.list")
  invoke("aggregate_pgen_susie.R", "missing metadata", c(aggregate_args, bad_args,
    "--output-dir", file.path(root, "missing metadata")), "readable")
  bad_args <- metadata_args
  bad_args[[4L]] <- write_list(output_triplets$pgen[-1L], "short pgens")
  invoke("aggregate_pgen_susie.R", "short metadata", c(aggregate_args, bad_args,
    "--output-dir", file.path(root, "short metadata")), "length")
  bad_args[[4L]] <- write_list(replace(output_triplets$pgen, 2L, "gs://bucket/other.pgen"), "conflicting pgens")
  invoke("aggregate_pgen_susie.R", "conflicting metadata", c(aggregate_args, bad_args,
    "--output-dir", file.path(root, "conflicting metadata")), "triplet")
  short_args <- unlist(lapply(names(output_metadata), function(name) {
    c(paste0("--", gsub("_", "-", name)), write_list(output_metadata[[name]][-1L], paste0("short ", name)))
  }), use.names = FALSE)
  # Output metadata can describe fewer entries or a different chromosome order
  # from the outcome index. It must not filter the supplied scientific tables.
  short_output <- file.path(root, "missing batch metadata")
  invoke("aggregate_pgen_susie.R", "missing batch metadata", c(aggregate_args, short_args,
    "--output-dir", short_output))
  short_manifest <- readPgenManifest(file.path(short_output, "pgen_manifest.tsv"))
  stopifnot(identical(short_manifest, tibble::tibble(chromosome = c("2", "3", "4", "5"),
    pgen = output_triplets$pgen[2:5], pvar = output_triplets$pvar[2:5],
    psam = output_triplets$psam[2:5])))
  bad_args <- metadata_args
  bad_args[[2L]] <- write_list(c("2", "3", "2", "4", "5"), "wrong batch chromosomes")
  reordered_output <- file.path(root, "wrong batch chromosomes")
  invoke("aggregate_pgen_susie.R", "wrong batch chromosomes", c(aggregate_args, bad_args,
    "--output-dir", reordered_output))
  reordered_manifest <- readPgenManifest(file.path(reordered_output, "pgen_manifest.tsv"))
  stopifnot(identical(reordered_manifest, tibble::tibble(chromosome = c("2", "3", "4", "5"),
    pgen = output_triplets$pgen[c(1L, 2L, 4L, 5L)],
    pvar = output_triplets$pvar[c(1L, 2L, 4L, 5L)],
    psam = output_triplets$psam[c(1L, 2L, 4L, 5L)])))
  for (output in c(short_output, reordered_output)) {
    for (name in c("susie.parquet", "credible_sets.parquet", "outcome_index.tsv", "covariate_explained_variants.tsv")) {
      stopifnot(identical(digest::digest(file = file.path(plain, name), algo = "sha256"),
                         digest::digest(file = file.path(output, name), algo = "sha256")))
    }
  }
  stopifnot(!file.exists("PWNED"))
  passed <- TRUE
  message("PGEN QC metadata tests passed: unique selected triplets, batch indices, CLI validation, and reusable URI manifests")
  invisible(TRUE)
}
