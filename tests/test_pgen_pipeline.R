# Real installed-code PLINK/PGEN pipeline tests. Run from the repository root.
benchmark_library <- Sys.getenv("BENCH_LIBRARY", unset = "")
if (nzchar(benchmark_library)) .libPaths(c(benchmark_library, .libPaths()))
args <- commandArgs(trailingOnly = TRUE)
repo_dir <- normalizePath(if (length(args)) args[[1L]] else getwd(), mustWork = TRUE)
functions_dir <- Sys.getenv("SUSIER_FUNCTIONS_PATH", file.path(repo_dir, "R", "utils"))
scripts_dir <- Sys.getenv("SUSIER_SCRIPTS_PATH", file.path(repo_dir, "R", "scripts"))
stopifnot("The batch CLI must be installed" =
            file.exists(file.path(scripts_dir, "run_pgen_susie_batch.R")))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(stringr))
source(file.path(scripts_dir, "run_pgen_susie_batch.R"))
source(file.path(repo_dir, "tests", "test_pgen_batch_inputs.R"))

run_pipeline_checks <- function() {
  previous_options <- options(warning.length = 8170)
  on.exit(options(previous_options), add = TRUE)
  plink <- Sys.getenv("PLINK2", unset = Sys.which("plink2"))
  stopifnot("Set PLINK2 to a real PLINK 2 executable" = nzchar(plink))
  directory <- tempfile("pgen pipeline 'quoted' $literal ", tmpdir = dirname(tempdir()))
  dir.create(directory)
  passed <- FALSE
  on.exit(if (passed) unlink(directory, recursive = TRUE) else
    message("Pipeline fixture and child logs retained at: ", directory), add = TRUE)
  rscript <- file.path(R.home("bin"), "Rscript")
  invoke <- function(script, arguments, label, expected_error = NULL) {
    log <- file.path(directory, paste0(label, ".log"))
    status <- suppressWarnings(system2(rscript, shQuote(c(file.path(scripts_dir, script), arguments)),
                                      stdout = log, stderr = log))
    output_lines <- readLines(log, warn = FALSE)
    output <- paste(output_lines, collapse = "\n")
    if (is.null(expected_error)) {
      if (!identical(status, 0L)) {
        cat(tail(output_lines, 80L), sep = "\n", file = stderr())
        stop(label, " failed; child log: ", log, call. = FALSE)
      }
    } else {
      stopifnot(status != 0L, grepl(expected_error, output, fixed = TRUE))
    }
    invisible(output)
  }
  set.seed(4821L)
  sample_ids <- as.character(1000L + rev(seq_len(120L)))
  positions <- seq.int(500L, 5500L, by = 5L)
  # Keep two distinct records at the same position; neither may be deduplicated.
  positions <- sort(c(positions, 1100L))
  refs <- rep("A", length(positions))
  alts <- rep("G", length(positions))
  alts[which(positions == 1100L)[[2L]]] <- "C"
  variant_ids <- paste0("chr1_", positions, "_", refs, "_", alts)
  hardcalls <- matrix(rbinom(length(positions) * length(sample_ids), 2L, 0.3),
                       nrow = length(positions), dimnames = list(variant_ids, sample_ids))
  hardcalls[sample(length(hardcalls), 80L)] <- NA_integer_
  hardcalls[positions == 1000L, ] <- 0L
  hardcalls[positions == 1000L, 1L] <- NA_integer_
  hardcalls[positions == 1005L, ] <- NA_integer_
  pgen <- file.path(directory, "input.pgen")
  pvar <- file.path(directory, "input.pvar")
  psam <- file.path(directory, "input.psam")
  writer <- pgenlibr::NewPgenWriter(pgen, length(sample_ids), nrow(hardcalls), nonref_flags = FALSE)
  for (i in seq_len(nrow(hardcalls))) pgenlibr::AppendBiallelic(writer, hardcalls[i, ])
  pgenlibr::ClosePgenWriter(writer)
  readr::write_tsv(tibble::tibble(`#CHROM` = "1", POS = positions,
                                 ID = paste0("rs", seq_along(positions)), REF = refs, ALT = alts), pvar)
  readr::write_tsv(tibble::tibble(`#IID` = sample_ids), psam)
  value_at <- function(position) {
    value <- hardcalls[which(positions == position)[[1L]], ]
    value[is.na(value)] <- mean(value, na.rm = TRUE)
    value
  }
  phenotype_ids <- c("outcomeC", "outcomeB2", "outcomeA", "outcomeB1")
  centers <- c(5000L, 3000L, 1100L, 1175L)
  values <- rbind(rnorm(120L), 1.8 * value_at(3000L) + rnorm(120L),
                   2 * value_at(1120L) + rnorm(120L),
                   1.5 * value_at(1175L) + rnorm(120L))
  colnames(values) <- sample_ids
  expression_samples <- sample(sample_ids)
  bed <- tibble::tibble(`#chr` = "1", start = centers, end = centers + 1L,
                         phenotype_id = phenotype_ids) |>
    dplyr::bind_cols(tibble::as_tibble(values[, expression_samples], .name_repair = "minimal"))
  phenotype_bed <- file.path(directory, "phenotypes.bed")
  readr::write_tsv(bed, phenotype_bed)
  qtl <- file.path(directory, "qtl.tsv")
  readr::write_tsv(tibble::tibble(phenotype_id = phenotype_ids, qval = .01), qtl)
  gene_map <- file.path(directory, "gene_map.tsv")
  readr::write_tsv(tibble::tibble(phenotype_id = phenotype_ids,
                                 gene_id = c("gene:C", "geneB", "gene A/one", "geneB")), gene_map)
  covariates <- file.path(directory, "covariates.tsv")
  cov <- rbind(value_at(1100L), rnorm(120L))
  colnames(cov) <- sample_ids
  readr::write_tsv(dplyr::bind_cols(tibble::tibble(ID = c("c1", "c2")),
                                    tibble::as_tibble(cov)), covariates)
  sample_list <- file.path(directory, "samples.tsv")
  readr::write_tsv(tibble::tibble(sample_id = sample(sample_ids)), sample_list)
  chromosomes <- file.path(directory, "chromosomes.txt")
  writeLines("chr1", chromosomes)
  pgen_manifest <- file.path(directory, "pgen_manifest.tsv")
  readr::write_tsv(tibble::tibble(chromosome = "chr1", pgen = pgen, pvar = pvar, psam = psam),
                    pgen_manifest)
  # The automatic cohort excludes both genotype-only and phenotype-only IDs.
  intersection_bed <- dplyr::select(bed, -dplyr::all_of(tail(sample_ids, 1L))) |>
    dplyr::mutate(trait_only = c(0, 1, 2, 3))
  intersection_bed_file <- file.path(directory, "intersection.bed")
  readr::write_tsv(intersection_bed, intersection_bed_file)
  intersection_ids <- sample_ids[sample_ids %in% names(intersection_bed)[-(1:4)]]
  stopifnot(length(intersection_ids) == 119L)
  intersection_sample_list <- file.path(directory, "intersection_samples.tsv")
  readr::write_tsv(tibble::tibble(sample_id = intersection_ids), intersection_sample_list)
  # Chromosome structs need not be in genomic order. The planner reports the
  # input array index, while grouping and the original region text stay stable.
  index_plan <- file.path(directory, "index_plan")
  planPgenBatches(phenotype_bed, qtl, c("chr2", "1"), 100, index_plan,
                  batch_size = 2L, gene_map = gene_map)
  stopifnot(identical(jsonlite::read_json(
    file.path(index_plan, "batch_chromosome_indices.json"), simplifyVector = TRUE),
    c(1L, 1L)))
  base_batch_args <- c("--pgen", pgen, "--pvar", pvar, "--psam", psam,
                       "--qtl-results", qtl, "--covariates", covariates,
                       "--cis-distance", "100",
                       "--base-seed", "82", "--plink2", plink)
  additional_bed <- file.path(directory, "additional.tsv")
  additional_ids <- c("chr1_3011_A_T", "chr2_1111_A_T", "chr1_1181_A_T")
  additional_samples <- rev(sample_ids[-1L])
  additional_values <- matrix(rbinom(3L * length(additional_samples), 2L, .3), nrow = 3L,
    dimnames = list(additional_ids, additional_samples))
  # The chr2 feature has the strong causal profile for the chr1 outcomeA. Its
  # duplicate cis profile keeps both records in the same high-purity set.
  additional_values[2L, ] <- value_at(1120L)[additional_samples]
  additional_values[1L, 3L] <- NA_real_
  additional_values[3L, 5L] <- -1
  readr::write_tsv(dplyr::bind_cols(tibble::tibble(
    variant_id = additional_ids, phenotype_id = c("outcomeB2", "outcomeA", "outcomeB1"),
    feature_type = "SV"), tibble::as_tibble(additional_values)), additional_bed, na = "NA")
  run_layout <- function(size, missing_code, additional = NULL, sample_mode = "provided") {
    label <- paste0("layout_", size, "_", gsub("-", "minus", missing_code, fixed = TRUE),
                    if (is.null(additional)) "" else "_additional",
                    if (sample_mode == "provided") "" else paste0("_", sample_mode))
    phenotype_input <- if (sample_mode == "provided") phenotype_bed else intersection_bed_file
    sample_args <- switch(sample_mode,
      provided = c("--sample-list", sample_list),
      inferred = character(),
      explicit_intersection = c("--sample-list", intersection_sample_list),
      stop("Unknown sample mode"))
    plan_dir <- file.path(directory, label, "plan")
    invoke("plan_pgen_susie.R", c(
      "--phenotype-bed", phenotype_input, "--qtl-results", qtl,
      "--pgen-manifest", pgen_manifest, "--cis-distance", "100",
      "--batch-size", as.character(size), "--gene-map", gene_map,
      "--output-dir", plan_dir), paste0(label, "_plan"))
    stopifnot("Planner must emit typed chromosome index metadata" =
                file.exists(file.path(plan_dir, "batch_chromosome_indices.json")))
    stopifnot(all(jsonlite::read_json(file.path(plan_dir, "batch_chromosome_indices.json"),
                                     simplifyVector = TRUE) == 0L))
    normalized <- readLines(file.path(plan_dir, "chromosome_pgens.tsv"), warn = FALSE)
    stopifnot(identical(normalized, paste("chr1", pgen, pvar, psam, sep = "\t")))
    worklists <- readLines(file.path(plan_dir, "batch_worklists.txt"))
    beds <- readLines(file.path(plan_dir, "batch_phenotype_beds.txt"))
    batch_dirs <- character(length(worklists))
    fits <- list()
    for (i in seq_along(worklists)) {
      batch_dirs[[i]] <- file.path(directory, label, paste0("run_", i))
      invoke("run_pgen_susie_batch.R", c(base_batch_args,
        sample_args,
        if (!is.null(additional)) c("--additional-genotypes", additional),
        "--worklist", worklists[[i]], "--phenotype-bed", beds[[i]],
        paste0("--missing-code=", missing_code), "--output-dir", batch_dirs[[i]]),
        paste0(label, "_batch_", i))
      work <- readr::read_tsv(worklists[[i]], show_col_types = FALSE)
      source_reader <- openPgenReader(pgen, pvar, psam)
      subset_reader <- openPgenReader(paste0(file.path(batch_dirs[[i]], "batch.genotypes"), ".pgen"),
                                       paste0(file.path(batch_dirs[[i]], "batch.genotypes"), ".pvar"),
                                       paste0(file.path(batch_dirs[[i]], "batch.genotypes"), ".psam"))
      keep <- Reduce(`|`, purrr::map(seq_len(nrow(work)), ~
        positions >= work$window_start[[.x]] & positions <= work$window_end[[.x]]))
      stopifnot(identical(subset_reader$variants$variant_id, variant_ids[keep]),
                identical(subset_reader$sample_ids, sample_ids))
      for (j in seq_len(nrow(work))) {
        actual <- readPgenWindow(subset_reader, "1", work$window_start[[j]], work$window_end[[j]])
        expected <- readPgenWindow(source_reader, "1", work$window_start[[j]], work$window_end[[j]])
        stopifnot(identical(actual, expected))
      }
      closePgenReader(source_reader)
      closePgenReader(subset_reader)
      rds_files <- list.files(file.path(batch_dirs[[i]], "genes"), pattern = "_susie[.]rds$",
                              full.names = TRUE, recursive = TRUE)
      for (rds_file in rds_files) {
        fitted <- readRDS(rds_file)
        fits[names(fitted)] <- fitted
        stopifnot(length(fitted) %in% c(1L, 2L))
      }
      gene_files <- list.files(file.path(batch_dirs[[i]], "genes"), recursive = TRUE)
      stopifnot(sum(grepl("[.]provenance[.]json$", gene_files)) == dplyr::n_distinct(work$gene_id))
    }
    lists <- c(summary = "batch.susie.parquet", cs = "batch.credible_sets.parquet",
                index = "batch.outcome_index.tsv", diagnostic = "batch.covariate_explained_variants.tsv")
    list_paths <- purrr::imap_chr(lists, function(name, key) {
      path <- file.path(directory, paste0(label, "_", key, ".files"))
      writeLines(file.path(batch_dirs, name), path)
      path
    })
    gene_files <- unlist(lapply(batch_dirs, function(path) {
      list.files(file.path(path, "genes"), recursive = TRUE, full.names = TRUE)
    }), use.names = FALSE)
    metadata_patterns <- c(rds = "_susie[.]rds$", parquet = "[.]parquet$",
      tsv = "[.]tsv$", log = "[.]log$", provenance = "[.]provenance[.]json$")
    metadata_args <- unlist(purrr::imap(metadata_patterns, function(pattern, role) {
      path <- file.path(directory, paste0(label, "_gene_", role, ".uris"))
      # Reversed file order must not change which gene owns each output.
      writeLines(rev(gene_files[grepl(pattern, gene_files)]), path)
      c(paste0("--gene-", role, "-uris-file"), path)
    }), use.names = FALSE)
    aggregate_dir <- file.path(directory, label, "aggregate")
    invoke("aggregate_pgen_susie.R", c(
      "--expected-outcomes", file.path(plan_dir, "expected_outcomes.tsv"),
      "--summary-files", list_paths[["summary"]], "--credible-set-files", list_paths[["cs"]],
      "--index-files", list_paths[["index"]], "--diagnostic-files", list_paths[["diagnostic"]],
      metadata_args, "--output-dir", aggregate_dir), paste0(label, "_aggregate"))
    manifests <- list.files(file.path(aggregate_dir, "gene_manifests"), full.names = TRUE)
    stopifnot(identical(basename(manifests), "chr1.gene_outputs.tsv"))
    manifest <- readr::read_tsv(manifests[[1L]],
      col_types = readr::cols(.default = readr::col_character()), show_col_types = FALSE)
    stopifnot(nrow(manifest) == 3L,
      setequal(manifest$gene_id, c("gene A/one", "geneB", "gene:C")),
      all(manifest$chromosome == "1"),
      setequal(manifest$rds_uri, gene_files[grepl("_susie[.]rds$", gene_files)]))
    manifest_uris <- unlist(manifest[grepl("_uri$", names(manifest))], use.names = FALSE)
    stopifnot(length(manifest_uris) == 27L, setequal(manifest_uris, gene_files),
      all(file.exists(manifest_uris)))
    index <- readr::read_tsv(file.path(aggregate_dir, "outcome_index.tsv"), show_col_types = FALSE)
    stopifnot(nrow(index) == 4L, setequal(index$outcome_id, phenotype_ids),
              !any(index$status == "FAILED"),
              identical(names(fits[sort(names(fits))]), sort(phenotype_ids)))
    provenance <- lapply(gene_files[grepl("[.]provenance[.]json$", gene_files)],
      jsonlite::read_json, simplifyVector = TRUE)
    expected_count <- if (sample_mode == "provided") 120L else 119L
    expected_mode <- if (sample_mode == "inferred") "genotype_phenotype_intersection" else "provided_list"
    for (record in provenance) {
      stopifnot(record$settings$sample_count == expected_count,
        identical(record$settings$sample_selection, expected_mode),
        identical(record$inputs$sample_list$present, sample_mode != "inferred"))
    }
    list(fits = fits[sort(names(fits))], index = index, provenance = provenance,
         summary = arrow::read_parquet(file.path(aggregate_dir, "susie.parquet")),
         credible_sets = arrow::read_parquet(file.path(aggregate_dir, "credible_sets.parquet")))
  }
  # These checks fail before a missing-file/URI can reach PLINK.
  invoke("run_pgen_susie_batch.R", c(base_batch_args, "--pgen", "gs://bucket/unlocalized.pgen",
    "--worklist", gene_map, "--phenotype-bed", phenotype_bed, "--missing-code=NA",
    "--output-dir", file.path(directory, "invalid")), "cloud_path", "localization error")
  invoke("run_pgen_susie_batch.R", c(base_batch_args, "--sample-list", "gs://bucket/unlocalized.samples",
    "--worklist", gene_map, "--phenotype-bed", phenotype_bed, "--missing-code=NA",
    "--output-dir", file.path(directory, "invalid_samples")), "cloud_samples", "localization error")
  additional_pair <- run_layout(2L, "NA", additional_bed)
  additional_single <- run_layout(1L, "NA", additional_bed)
  stopifnot("Additional features must not depend on batch size" =
    identical(additional_pair$fits, additional_single$fits))
  stopifnot("The gene-owned chr2 feature must reach the aggregated credible-set summary" =
    any(additional_pair$summary$variant == "chr2_1111_A_T" &
          additional_pair$summary$outcome_id == "outcomeA" &
          additional_pair$summary$chromosome == "2"))
  for (id in phenotype_ids) {
    expected_additional <- if (id == "outcomeA") additional_ids[2L] else
      if (id %in% c("outcomeB1", "outcomeB2")) additional_ids[c(1L, 3L)] else character()
    actual_ids <- additional_pair$fits[[id]]$variant_id
    stopifnot(setequal(intersect(actual_ids, additional_ids), expected_additional),
      !any(setdiff(additional_ids, expected_additional) %in% actual_ids))
  }
  first <- run_layout(2L, "NA")
  single <- run_layout(1L, "NA")
  hundred <- run_layout(100L, "NA")
  stopifnot(identical(first$fits, single$fits), identical(first$fits, hundred$fits))
  first_minus <- run_layout(2L, "-1")
  single_minus <- run_layout(1L, "-1")
  stopifnot(identical(first_minus$fits, single_minus$fits),
    "Native PGEN fits must not depend on the CLI missing-code flag" =
      identical(first$fits, first_minus$fits))
  compare_summary <- function(left, right) {
    columns <- setdiff(names(left), "batch_id")
    left <- dplyr::arrange(left[columns], gene_id, outcome_id, dplyr::across(dplyr::everything()))
    right <- dplyr::arrange(right[columns], gene_id, outcome_id, dplyr::across(dplyr::everything()))
    stopifnot(identical(as.data.frame(left), as.data.frame(right)))
  }
  compare_summary(first$summary, single$summary)
  compare_summary(first$credible_sets, single$credible_sets)
  for (additional in list(NULL, additional_bed)) {
    inferred <- run_layout(100L, "NA", additional, sample_mode = "inferred")
    explicit <- run_layout(100L, "NA", additional, sample_mode = "explicit_intersection")
    stopifnot("Automatic and explicit intersection cohorts must give identical fits" =
      identical(inferred$fits, explicit$fits))
    compare_summary(inferred$summary, explicit$summary)
    compare_summary(inferred$credible_sets, explicit$credible_sets)
    inferred_hashes <- vapply(inferred$provenance, function(record) record$settings$sample_ids_sha256,
      character(1))
    explicit_hashes <- vapply(explicit$provenance, function(record) record$settings$sample_ids_sha256,
      character(1))
    stopifnot(identical(inferred_hashes, explicit_hashes))
  }
  # Native PGEN-only reads use NA missing semantics for either CLI flag.
  # A -1 Tabix reference would keep all-missing and zero-ALT rows during the
  # pre-imputation filter. Use NA to compare the same retained variants and seeds.
  for (missing_code in c("NA", "-1")) {
    doses <- hardcalls
    dosage_tsv <- file.path(directory, paste0("reference_", missing_code, ".tsv"))
    readr::write_tsv(dplyr::bind_cols(
      tibble::tibble(CHROM = "1", POS = positions, REF = refs, ALT = alts),
      tibble::as_tibble(doses)), dosage_tsv, na = "NA")
    dosage_bgz <- paste0(dosage_tsv, ".bgz")
    Rsamtools::bgzip(dosage_tsv, dosage_bgz, overwrite = TRUE)
    Rsamtools::indexTabix(dosage_bgz, seq = 1L, start = 2L, end = 2L, skip = 1L)
    data <- LoadData(list(expression_matrix = phenotype_bed, phenotype_list = qtl,
      covariates = covariates, genotype_matrix = dosage_bgz, sample_meta = sample_list,
      cisdistance = 100, out_prefix = "reference", MAF = 0))
    se <- makePgenBatchExperiment(data)
    legacy_se <- eQTLUtils::makeSummarizedExperimentFromCountMatrix(
      dplyr::select(data$expression_matrix, -1, -2, -3), data$phenotype_meta,
      data$sample_metadata, quant_method = "gene_counts", reformat = FALSE)
    stopifnot(identical(se, legacy_se))
    expected_fits <- if (missing_code == "NA") first$fits else first_minus$fits
    for (id in phenotype_ids) {
      gene_id <- readr::read_tsv(gene_map, show_col_types = FALSE)$gene_id[match(id, phenotype_ids)]
      reference <- .pgenGeneWithSeed(geneOutcomeSeed(82L, gene_id, id),
        finemapPhenotype(id, se, dosage_bgz, data$covariates_matrix, 100))
      stopifnot(identical(reference, expected_fits[[id]]))
      if (missing_code == "NA") {
        # Recreate the old per-gene task: its BED contains every selected outcome
        # for this gene, so LoadData supplies their union to each outcome.
        gene_ids <- phenotype_ids[readr::read_tsv(gene_map, show_col_types = FALSE)$gene_id == gene_id]
        gene_bed <- file.path(directory, "legacy_gene.bed")
        readr::write_tsv(dplyr::filter(bed, phenotype_id %in% gene_ids), gene_bed)
        gene_data <- LoadData(list(expression_matrix = gene_bed, phenotype_list = qtl,
          covariates = covariates, genotype_matrix = dosage_bgz, sample_meta = sample_list,
          cisdistance = 100, out_prefix = "reference", MAF = 0,
          AdditionalGenotypesBed = additional_bed))
        gene_se <- makePgenBatchExperiment(gene_data)
        gene_reference <- .pgenGeneWithSeed(geneOutcomeSeed(82L, gene_id, id),
          finemapPhenotype(id, gene_se, dosage_bgz, gene_data$covariates_matrix, 100,
            additional_genotypes = gene_data$additional_genotypes))
        stopifnot("Gene-local additional features must match the legacy Tabix fit" =
          identical(gene_reference, additional_pair$fits[[id]]))
      }
    }
  }
  message("PGEN pipeline tests passed: real PLINK unions, exact Tabix fit parity, ",
          "batch sizes 1/2/100, both missing codes, separate gene artifacts, complete aggregation, ",
          "chromosome gene manifests, PGEN input manifest, and inferred sample intersections")
  passed <- TRUE
}
run_pipeline_checks()
