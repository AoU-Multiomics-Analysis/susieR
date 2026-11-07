# Small fixtures check gene identity and preserve artifact locations as metadata.
run_gene_manifest_tests <- function(repo_dir) {
  suppressPackageStartupMessages(library(dplyr))
  functions_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH",
                               unset = file.path(repo_dir, "R", "utils"))
  source(file.path(functions_path, "PgenBatchFunctions.R"))
  stopifnot("The chromosome gene manifest builder must be available" =
              exists("writeChromosomeGeneManifests", mode = "function"))
  root <- tempfile("gene manifests 'quoted' ")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  index_file <- file.path(root, "combined outcome index.tsv")
  genes <- tibble::tibble(
    gene_id = c("gene A/α", "geneB", "gene:C"),
    gene_key = c(
      "gene-fc2a39a48c4febe822dc2d08b8a3ddf937ca9b39638fe0a3eef66f7d7959d3e5",
      "gene-7fb10082205fbb73dc840c440f2ce171fa5781f47e6e6481c5e05942ba46c0a2",
      "gene-49330230f76199b695bbbdef19f5e9334b08e958e8bc92c08a310187bfeaacea"))
  index <- tibble::tibble(
    batch_id = c("batch_2", "batch_1", "batch_2", "batch_3"),
    chromosome = c("chr2", "1", "2", "2"),
    gene_id = genes$gene_id[c(1, 2, 1, 3)],
    gene_key = genes$gene_key[c(1, 2, 1, 3)],
    outcome_id = c("a1", "b1", "a2", "c1"),
    phenotype_id = c("trait1", "trait2", "trait3", "trait4"),
    seed = 1:4, status = c("CONVERGED", "NONCONVERGED", "NONCONVERGED", "CONVERGED"),
    converged = c(TRUE, FALSE, FALSE, TRUE), credible_set_count = c(1L, 0L, 0L, 0L),
    reason = c("", "SuSiE did not converge", "SuSiE did not converge", ""))
  readr::write_tsv(index, index_file)
  suffixes <- c(
    rds_uri = "_susie.rds", summary_parquet_uri = ".parquet",
    lbf_parquet_uri = ".lbf_variable.parquet", full_parquet_uri = ".full_susie.parquet",
    covariate_explained_variants_uri = ".covariate_explained_variants.tsv",
    variant_position_summary_uri = ".variant_position_summary.tsv",
    outcome_index_uri = ".outcome_index.tsv", log_uri = ".log",
    provenance_uri = ".provenance.json")
  # These locations do not exist. The builder must never open them.
  prefixes <- paste0(c("gs://missing-bucket/workflow path/", "/missing local output/",
                       "gs://missing-bucket/another batch/"), genes$gene_key, "/", genes$gene_key)
  artifacts <- lapply(suffixes, function(suffix) paste0(prefixes, suffix))
  arguments <- list(
    outcome_index_file = index_file,
    gene_rds_uris = rev(artifacts$rds_uri),
    gene_parquet_uris = rev(unlist(artifacts[c("full_parquet_uri", "summary_parquet_uri",
                                               "lbf_parquet_uri")], use.names = FALSE)),
    gene_tsv_uris = rev(unlist(artifacts[c("outcome_index_uri", "covariate_explained_variants_uri",
                                           "variant_position_summary_uri")], use.names = FALSE)),
    gene_log_uris = artifacts$log_uri[c(2, 3, 1)],
    gene_provenance_uris = artifacts$provenance_uri[c(3, 1, 2)],
    output_dir = file.path(root, "output manifests"))
  files <- do.call(writeChromosomeGeneManifests, arguments)
  stopifnot(is.character(files), identical(names(files), c("1", "2")),
            identical(basename(unname(files)), c("chr1.gene_outputs.tsv", "chr2.gene_outputs.tsv")),
            all(file.exists(files)))
  read_manifest <- function(path) readr::read_tsv(path,
    col_types = readr::cols(.default = readr::col_character()), show_col_types = FALSE)
  chr1 <- read_manifest(files[["1"]])
  chr2 <- read_manifest(files[["2"]])
  stopifnot(identical(names(chr1), c("chromosome", "gene_id", "gene_key", "batch_id", names(suffixes))),
            identical(chr1$gene_id, "geneB"), identical(chr2$gene_id, c("gene A/α", "gene:C")),
            identical(chr1$chromosome, "1"), identical(chr2$chromosome, c("2", "2")),
            identical(chr2$batch_id, c("batch_2", "batch_3")))
  combined <- dplyr::bind_rows(chr2[1, ], chr1, chr2[2, ])
  stopifnot(identical(combined$gene_id, genes$gene_id),
            identical(combined$gene_key, genes$gene_key))
  for (column in names(artifacts)) {
    stopifnot("Every artifact URI must remain unchanged and belong to its gene" =
                identical(combined[[column]], artifacts[[column]]))
  }
  # Gene IDs and output-only URI text may contain literal quotes.
  quoted_gene <- 'gene"quoted'
  quoted_key <- paste0("gene-", digest::digest(charToRaw(quoted_gene),
                                               algo = "sha256", serialize = FALSE))
  quoted_index <- index %>% mutate(
    gene_id = replace(gene_id, gene_id == genes$gene_id[[1L]], quoted_gene),
    gene_key = replace(gene_key, gene_key == genes$gene_key[[1L]], quoted_key))
  readr::write_tsv(quoted_index, index_file, quote = "needed")
  quoted_arguments <- arguments
  quoted_arguments$output_dir <- file.path(root, "quoted manifests")
  for (group in c("gene_rds_uris", "gene_parquet_uris", "gene_tsv_uris",
                  "gene_log_uris", "gene_provenance_uris")) {
    quoted_arguments[[group]] <- gsub(genes$gene_key[[1L]], quoted_key,
                                      quoted_arguments[[group]], fixed = TRUE)
    quoted_arguments[[group]] <- gsub("workflow path", 'workflow "quoted" path',
                                      quoted_arguments[[group]], fixed = TRUE)
  }
  quoted_files <- do.call(writeChromosomeGeneManifests, quoted_arguments)
  quoted_chr2 <- read_manifest(quoted_files[["2"]])
  stopifnot(identical(quoted_chr2$gene_id[[1L]], quoted_gene),
            grepl('workflow "quoted" path', quoted_chr2$rds_uri[[1L]], fixed = TRUE))
  readr::write_tsv(index, index_file)
  # Index reads and output writes are allowed. Artifact filesystem probes are not.
  guarded_builder <- writeChromosomeGeneManifests
  guard <- new.env(parent = environment(guarded_builder))
  for (helper in c("file.exists", "file.access", "file.info", "dir.exists", "normalizePath",
                   "readLines", "readRDS", "readBin", "file", "gzfile")) {
    assign(helper, function(...) stop("Artifact metadata reached a filesystem helper", call. = FALSE),
           envir = guard)
  }
  environment(guarded_builder) <- guard
  guarded_arguments <- arguments
  guarded_arguments$output_dir <- file.path(root, "guarded manifests")
  guarded_files <- do.call(guarded_builder, guarded_arguments)
  for (chromosome in names(files)) {
    stopifnot(identical(readLines(guarded_files[[chromosome]]), readLines(files[[chromosome]])))
  }
  # Rejected inputs must not produce a partial set of chromosome manifests.
  check_rejected <- function(candidate, pattern) {
    candidate$output_dir <- file.path(root, "rejected manifests")
    unlink(candidate$output_dir, recursive = TRUE)
    result <- tryCatch(do.call(writeChromosomeGeneManifests, candidate), error = identity)
    if (!inherits(result, "error")) return("accepted invalid input")
    if (!grepl(pattern, conditionMessage(result), ignore.case = TRUE)) {
      return(paste("unexpected error:", conditionMessage(result)))
    }
    if (length(list.files(candidate$output_dir, pattern = "gene_outputs[.]tsv$"))) {
      return("wrote partial manifests before rejecting input")
    }
    NULL
  }
  invalid_indices <- list(
    empty_index = index[FALSE, ],
    missing_identity_column = dplyr::select(index, -gene_key),
    duplicate_outcome_row = dplyr::bind_rows(index, index[1, ]),
    duplicate_outcome_id = dplyr::mutate(index, outcome_id = c("a1", "a1", "a2", "c1")),
    duplicate_phenotype_id = dplyr::mutate(index, phenotype_id = c("trait1", "trait1", "trait3", "trait4")),
    absent_gene_id = dplyr::mutate(index, gene_id = replace(gene_id, 1, NA_character_)),
    blank_batch_id = dplyr::mutate(index, batch_id = replace(batch_id, 1, "")),
    multiline_gene_id = dplyr::mutate(index, gene_id = replace(gene_id, 1, "bad\nID")),
    wrong_gene_key = dplyr::mutate(index, gene_key = replace(gene_key, 1, genes$gene_key[[2]])),
    split_gene_batch = dplyr::mutate(index, batch_id = replace(batch_id, 3, "batch_4")),
    split_gene_chromosome = dplyr::mutate(index, chromosome = replace(chromosome, 3, "chr1")),
    unsupported_chromosome = dplyr::mutate(index, chromosome = replace(chromosome, 1, "chr23")),
    failed_outcome = dplyr::mutate(index, status = replace(status, 1, "FAILED")),
    unknown_status = dplyr::mutate(index, status = replace(status, 1, "SKIPPED")),
    missing_status = dplyr::mutate(index, status = replace(status, 1, NA_character_)),
    inconsistent_convergence = dplyr::mutate(index, converged = replace(converged, 1, FALSE)),
    missing_convergence = dplyr::mutate(index, converged = replace(converged, 1, NA)))
  failures <- character()
  for (name in names(invalid_indices)) {
    candidate <- arguments
    readr::write_tsv(invalid_indices[[name]], index_file)
    failure <- check_rejected(candidate, "index|schema|column|duplicate|gene|batch|chromosome|autosome|status|converged|boolean")
    if (!is.null(failure)) failures[[name]] <- failure
  }
  readr::write_tsv(index, index_file)
  if (length(failures)) {
    stop("Invalid outcome indices were not rejected safely: ",
         paste(paste(names(failures), failures, sep = " = "), collapse = "; "), call. = FALSE)
  }
  invalid_artifacts <- list()
  changed <- arguments
  changed$gene_rds_uris <- changed$gene_rds_uris[-1]
  invalid_artifacts$missing_rds <- changed
  changed <- arguments
  changed$gene_parquet_uris <- changed$gene_parquet_uris[-1]
  invalid_artifacts$missing_parquet_role <- changed
  changed <- arguments
  changed$gene_rds_uris <- c(changed$gene_rds_uris, changed$gene_rds_uris[[1]])
  invalid_artifacts$duplicate_uri <- changed
  changed <- arguments
  changed$gene_rds_uris <- c(changed$gene_rds_uris, sub("missing-bucket", "another-bucket", changed$gene_rds_uris[[1]]))
  invalid_artifacts$duplicate_role_at_another_location <- changed
  changed <- arguments
  changed$gene_parquet_uris[[1]] <- "gs://missing-bucket/gene-foreign.parquet"
  invalid_artifacts$foreign_gene <- changed
  changed <- arguments
  changed$gene_tsv_uris[[1]] <- "gs://missing-bucket/unrelated.tsv"
  invalid_artifacts$unknown_suffix <- changed
  changed <- arguments
  changed$gene_log_uris[[1]] <- changed$gene_provenance_uris[[1]]
  changed$gene_provenance_uris[[1]] <- arguments$gene_log_uris[[1]]
  invalid_artifacts$wrong_artifact_group <- changed
  for (value in list(NA_character_, "", "bad\nlocation", "bad\tlocation", 123, NULL, character())) {
    changed <- arguments
    changed["gene_rds_uris"] <- list(value)
    invalid_artifacts[[paste0("malformed_rds_vector_", length(invalid_artifacts) + 1L)]] <- changed
  }
  failures <- character()
  for (name in names(invalid_artifacts)) {
    failure <- check_rejected(invalid_artifacts[[name]], "artifact|URI|metadata|role|group|duplicate|missing|foreign")
    if (!is.null(failure)) failures[[name]] <- failure
  }
  if (length(failures)) {
    stop("Invalid artifact metadata was not rejected safely: ",
         paste(paste(names(failures), failures, sep = " = "), collapse = "; "), call. = FALSE)
  }
  message("Chromosome gene manifest identity and URI preservation checks passed.")
}
