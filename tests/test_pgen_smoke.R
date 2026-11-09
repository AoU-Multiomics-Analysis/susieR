# One real installed-code planner, PLINK extraction, fit and aggregate check.
benchmark_library <- Sys.getenv("BENCH_LIBRARY", unset = "")
if (nzchar(benchmark_library)) .libPaths(c(benchmark_library, .libPaths()))
functions_dir <- Sys.getenv("SUSIER_FUNCTIONS_PATH", unset = "/opt/r/lib")
scripts_dir <- Sys.getenv("SUSIER_SCRIPTS_PATH", unset = "/opt/r/scripts")
source(file.path(functions_dir, "PgenFunctions.R"))

run_pgen_smoke <- function() {
  plink <- Sys.getenv("PLINK2", unset = Sys.which("plink2"))
  stopifnot("A real PLINK 2 executable is required" = nzchar(plink))
  scripts <- c("plan_pgen_susie.R", "filter_pgen_susie.R", "run_pgen_susie_batch.R", "aggregate_pgen_susie.R")
  stopifnot("The PGEN CLIs must be installed" =
              all(file.exists(file.path(scripts_dir, scripts))))
  directory <- tempfile("pgen smoke 'quoted' $literal ", tmpdir = dirname(tempdir()))
  dir.create(directory)
  passed <- FALSE
  on.exit(if (passed) unlink(directory, recursive = TRUE) else
    message("Smoke fixture and child logs retained at: ", directory), add = TRUE)
  rscript <- file.path(R.home("bin"), "Rscript")
  invoke <- function(script, arguments, label, expected_error = NULL) {
    message("CHECK ", label)
    log <- file.path(directory, paste0(label, ".log"))
    status <- suppressWarnings(system2(
      rscript, shQuote(c(file.path(scripts_dir, script), arguments)),
      stdout = log, stderr = log
    ))
    output <- readLines(log, warn = FALSE)
    if (is.null(expected_error)) {
      if (!identical(status, 0L)) {
        cat(tail(output, 80L), sep = "\n", file = stderr())
        stop(label, " failed; child log: ", log, call. = FALSE)
      }
    } else {
      stopifnot(status != 0L,
                any(grepl(expected_error, output, fixed = TRUE)))
    }
    invisible(output)
  }

  message("PREPARE 48 samples, 51 variants and one phenotype")
  set.seed(4821L)
  sample_ids <- as.character(1000L + rev(seq_len(48L)))
  positions <- seq.int(875L, 1125L, by = 5L)
  variant_ids <- paste0("chr1_", positions, "_A_G")
  hardcalls <- matrix(rbinom(length(positions) * length(sample_ids), 2L, .3),
    nrow = length(positions), dimnames = list(variant_ids, sample_ids))
  storage.mode(hardcalls) <- "integer"
  hardcalls[sample(length(hardcalls), 8L)] <- NA_integer_
  hardcalls[positions == 905L, 1:3] <- NA_integer_
  pgen <- file.path(directory, "input.pgen")
  pvar <- file.path(directory, "input.pvar")
  psam <- file.path(directory, "input.psam")
  writer <- pgenlibr::NewPgenWriter(pgen, length(sample_ids), nrow(hardcalls),
                                   nonref_flags = FALSE)
  for (i in seq_len(nrow(hardcalls))) pgenlibr::AppendBiallelic(writer, hardcalls[i, ])
  pgenlibr::ClosePgenWriter(writer)
  readr::write_tsv(tibble::tibble(`#CHROM` = "1", POS = positions,
    ID = paste0("rs", seq_along(positions)), REF = "A", ALT = "G"), pvar)
  readr::write_tsv(tibble::tibble(`#IID` = sample_ids), psam)

  phenotype_id <- "outcomeA"
  gene_id <- "gene A/one"
  causal <- hardcalls[positions == 1000L, ]
  causal[is.na(causal)] <- mean(causal, na.rm = TRUE)
  phenotype <- 2.5 * causal + rnorm(length(sample_ids))
  phenotype_samples <- sample(sample_ids)
  phenotype_bed <- file.path(directory, "phenotypes.bed")
  # Keep the quoted, padded raw ID. The planner must select its decoded ID.
  bed_lines <- c(
    paste(c("#chr", "start", "end", "phenotype_id", phenotype_samples), collapse = "\t"),
    paste(c("chr1", "1000", "1001", '  "outcomeA"  ',
            phenotype[match(phenotype_samples, sample_ids)]), collapse = "\t")
  )
  writeLines(bed_lines, phenotype_bed)
  qtl <- file.path(directory, "qtl.tsv")
  readr::write_tsv(tibble::tibble(phenotype_id = phenotype_id, qval = .01), qtl)
  gene_map <- file.path(directory, "gene_map.tsv")
  readr::write_tsv(tibble::tibble(phenotype_id = phenotype_id, gene_id = gene_id), gene_map)
  covariates <- file.path(directory, "covariates.tsv")
  cov <- matrix(rnorm(2L * length(sample_ids)), nrow = 2L,
                 dimnames = list(NULL, sample_ids))
  readr::write_tsv(dplyr::bind_cols(tibble::tibble(ID = c("c1", "c2")),
                                    tibble::as_tibble(cov)), covariates)
  pgen_manifest <- file.path(directory, "pgen_manifest.tsv")
  readr::write_tsv(tibble::tibble(chromosome = "chr1", pgen = pgen,
                                  pvar = pvar, psam = psam), pgen_manifest)

  plan_dir <- file.path(directory, "plan")
  invoke("plan_pgen_susie.R", c(
    "--phenotype-bed", phenotype_bed, "--qtl-results", qtl,
    "--pgen-manifest", pgen_manifest, "--cis-distance", "100",
    "--batch-size", "1", "--gene-map", gene_map, "--output-dir", plan_dir
  ), "planner")
  stopifnot("Planner must preserve the typed chromosome input index" =
    identical(jsonlite::read_json(file.path(plan_dir, "batch_chromosome_indices.json"),
                                  simplifyVector = TRUE), 0L))
  stopifnot(identical(readLines(file.path(plan_dir, "chromosome_pgens.tsv"), warn = FALSE),
                       paste("chr1", pgen, pvar, psam, sep = "\t")))
  worklists <- readLines(file.path(plan_dir, "batch_worklists.txt"), warn = FALSE)
  beds <- readLines(file.path(plan_dir, "batch_phenotype_beds.txt"), warn = FALSE)
  stopifnot(length(worklists) == 1L, length(beds) == 1L,
            identical(readLines(beds[[1L]], warn = FALSE), bed_lines))
  work <- readr::read_tsv(worklists[[1L]], show_col_types = FALSE)
  stopifnot(nrow(work) == 1L, identical(work$phenotype_id, phenotype_id),
            identical(work$gene_id, gene_id), work$window_start == 900L,
            work$window_end == 1100L)

  qc_dir <- file.path(directory, "qc")
  invoke("filter_pgen_susie.R", c("--pgen", pgen, "--pvar", pvar, "--psam", psam,
    "--mac", "0", "--min-call-rate", "0.95", "--plink2", plink,
    "--output-dir", qc_dir), "shared_qc")
  prepared <- paste0(file.path(qc_dir, "filtered"), c(".pgen", ".pvar", ".psam"))
  qc_report <- file.path(qc_dir, "qc.json")
  qc <- jsonlite::read_json(qc_report, simplifyVector = TRUE)
  stopifnot(qc$sample_count == 48L, qc$variant_count_before == 51L,
    qc$variant_count_after == 50L, qc$thresholds$min_call_rate == .95,
    identical(jsonlite::read_json(file.path(plan_dir, "batch_pgen_indices.json"),
                                  simplifyVector = TRUE), 0L))

  batch_dir <- file.path(directory, "batch")
  batch_args <- c(
    "--pgen", prepared[[1L]], "--pvar", prepared[[2L]], "--psam", prepared[[3L]],
    "--genotype-qc-report", qc_report,
    "--worklist", worklists[[1L]], "--phenotype-bed", beds[[1L]],
    "--qtl-results", qtl, "--covariates", covariates, "--cis-distance", "100",
    "--base-seed", "82", "--plink2", plink, "--missing-code=NA",
    "--output-dir", batch_dir
  )
  invoke("run_pgen_susie_batch.R", c(batch_args, "--pgen", "gs://bucket/input.pgen"),
         "unlocalized_input", "localization error")
  # No sample-list argument: use the genotype/phenotype intersection in PSAM order.
  invoke("run_pgen_susie_batch.R", batch_args, "batch")
  extracted_files <- paste0(file.path(batch_dir, "batch.genotypes"),
                             c(".pgen", ".pvar", ".psam"))
  stopifnot(all(file.exists(extracted_files)))
  reader <- openPgenReader(extracted_files[[1L]], extracted_files[[2L]], extracted_files[[3L]])
  on.exit(closePgenReader(reader), add = TRUE)
  keep <- positions >= 900L & positions <= 1100L & positions != 905L
  stopifnot(sum(keep) == 40L, sum(!keep) == 11L,
            identical(reader$variants$variant_id, variant_ids[keep]),
            identical(reader$sample_ids, sample_ids),
            identical(readPgenWindow(reader, "1", 900L, 1100L), hardcalls[keep, , drop = FALSE]))
  closePgenReader(reader)

  gene_files <- list.files(file.path(batch_dir, "genes"), recursive = TRUE, full.names = TRUE)
  rds_files <- gene_files[grepl("_susie[.]rds$", gene_files)]
  stopifnot(length(rds_files) == 1L)
  fits <- readRDS(rds_files[[1L]])
  stopifnot(identical(names(fits), phenotype_id))
  fitted <- fits[[1L]]
  stopifnot("The smoke test must contain a real fitted outcome" =
    inherits(fitted, "susie"), isTRUE(fitted$converged),
    identical(fitted$variant_id, variant_ids[keep]),
    length(fitted$pip) == sum(keep), all(is.finite(fitted$pip)),
    all(fitted$pip >= 0 & fitted$pip <= 1))
  provenance_files <- gene_files[grepl("[.]provenance[.]json$", gene_files)]
  stopifnot(length(provenance_files) == 1L)
  provenance <- jsonlite::read_json(provenance_files[[1L]], simplifyVector = TRUE)
  stopifnot(provenance$settings$sample_count == 48L,
    provenance$settings$min_call_rate == .95,
    provenance$settings$mac == 0,
    identical(provenance$genotype_qc$count_source, "pgen_hardcalls"),
    identical(provenance$settings$sample_selection, "genotype_phenotype_intersection"),
    identical(provenance$settings$sample_order, "psam"),
    identical(provenance$inputs$sample_list$present, FALSE),
    identical(provenance$settings$sample_ids_sha256,
      digest::digest(sample_ids, algo = "sha256", serializeVersion = 2)))

  batch_outputs <- c(summary = "batch.susie.parquet", cs = "batch.credible_sets.parquet",
    index = "batch.outcome_index.tsv", diagnostic = "batch.covariate_explained_variants.tsv")
  list_paths <- purrr::imap_chr(batch_outputs, function(filename, role) {
    path <- file.path(directory, paste0(role, ".files"))
    writeLines(file.path(batch_dir, filename), path)
    path
  })
  metadata_patterns <- c(rds = "_susie[.]rds$", parquet = "[.]parquet$",
    tsv = "[.]tsv$", log = "[.]log$", provenance = "[.]provenance[.]json$")
  metadata_args <- unlist(purrr::imap(metadata_patterns, function(pattern, role) {
    path <- file.path(directory, paste0("gene_", role, ".uris"))
    writeLines(gene_files[grepl(pattern, gene_files)], path)
    c(paste0("--gene-", role, "-uris-file"), path)
  }), use.names = FALSE)
  aggregate_dir <- file.path(directory, "aggregate")
  genotype_metadata <- c(pgen_chromosomes = "1", pgen_uris = prepared[[1L]],
                         pvar_uris = prepared[[2L]], psam_uris = prepared[[3L]])
  genotype_args <- unlist(purrr::imap(genotype_metadata, function(value, role) {
    path <- file.path(directory, paste0(role, ".metadata"))
    writeLines(value, path)
    c(paste0("--", gsub("_", "-", role), "-file"), path)
  }), use.names = FALSE)
  invoke("aggregate_pgen_susie.R", c(
    "--expected-outcomes", file.path(plan_dir, "expected_outcomes.tsv"),
    "--summary-files", list_paths[["summary"]], "--credible-set-files", list_paths[["cs"]],
    "--index-files", list_paths[["index"]], "--diagnostic-files", list_paths[["diagnostic"]],
    metadata_args, genotype_args, "--output-dir", aggregate_dir
  ), "aggregate")
  stopifnot(all(file.exists(file.path(aggregate_dir, c("susie.parquet", "credible_sets.parquet",
    "outcome_index.tsv", "covariate_explained_variants.tsv")))))
  index <- readr::read_tsv(file.path(aggregate_dir, "outcome_index.tsv"), show_col_types = FALSE)
  stopifnot(nrow(index) == 1L, identical(index$outcome_id, phenotype_id),
            identical(index$gene_id, gene_id), identical(index$status, "CONVERGED"),
            isTRUE(index$converged))
  summary <- arrow::read_parquet(file.path(aggregate_dir, "susie.parquet"))
  credible_sets <- arrow::read_parquet(file.path(aggregate_dir, "credible_sets.parquet"))
  stopifnot(nrow(summary) > 0L, nrow(credible_sets) > 0L,
            all(summary$outcome_id == phenotype_id),
            all(credible_sets$outcome_id == phenotype_id))
  manifests <- list.files(file.path(aggregate_dir, "gene_manifests"), full.names = TRUE)
  stopifnot(identical(basename(manifests), "chr1.gene_outputs.tsv"))
  manifest <- readr::read_tsv(manifests[[1L]],
    col_types = readr::cols(.default = readr::col_character()), show_col_types = FALSE)
  manifest_uris <- unlist(manifest[grepl("_uri$", names(manifest))], use.names = FALSE)
  stopifnot(nrow(manifest) == 1L, identical(manifest$gene_id, gene_id),
            identical(manifest$chromosome, "1"), length(manifest_uris) == 9L,
            setequal(manifest_uris, gene_files), all(file.exists(manifest_uris)))
  reusable <- readr::read_tsv(file.path(aggregate_dir, "pgen_manifest.tsv"),
    col_types = readr::cols(.default = readr::col_character()), show_col_types = FALSE)
  stopifnot(nrow(reusable) == 1L, identical(reusable$chromosome, "1"),
    identical(unlist(reusable[c("pgen", "pvar", "psam")], use.names = FALSE), prepared))
  passed <- TRUE
  message("PGEN smoke passed: one planner, PLINK extraction, fitted outcome and aggregate")
}

run_pgen_smoke()
