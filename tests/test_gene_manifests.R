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
    col_types = readr::cols(.default = readr::col_character()), na = "NA", show_col_types = FALSE)
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
  # Output metadata remains available for repeated and failed outcome rows.
  repeated_index <- bind_rows(index, index[1L, ]) %>% mutate(status = "FAILED", converged = FALSE)
  readr::write_tsv(repeated_index, index_file)
  repeated_arguments <- arguments
  repeated_arguments$output_dir <- file.path(root, "repeated failed manifests")
  repeated_files <- do.call(writeChromosomeGeneManifests, repeated_arguments)
  repeated_rows <- bind_rows(lapply(repeated_files, read_manifest))
  stopifnot(nrow(repeated_rows) == 3L, setequal(repeated_rows$gene_id, genes$gene_id))
  readr::write_tsv(index, index_file)
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
  # Gene metadata is written without checks on scientific status or identity.
  accepted_indices <- list(
    empty_index = index[FALSE, ],
    missing_identity_column = dplyr::select(index, -gene_key),
    duplicate_outcome_row = dplyr::bind_rows(index, index[1, ]),
    duplicate_outcome_id = dplyr::mutate(index, outcome_id = c("a1", "a1", "a2", "c1")),
    duplicate_phenotype_id = dplyr::mutate(index, phenotype_id = c("trait1", "trait1", "trait3", "trait4")),
    absent_gene_id = dplyr::mutate(index, gene_id = replace(gene_id, 1, NA_character_)),
    blank_batch_id = dplyr::mutate(index, batch_id = replace(batch_id, 1, "")),
    multiline_gene_id = dplyr::mutate(index, gene_id = replace(gene_id, 1, "gene\nID")),
    wrong_gene_key = dplyr::mutate(index, gene_key = replace(gene_key, 1, genes$gene_key[[2]])),
    split_gene_batch = dplyr::mutate(index, batch_id = replace(batch_id, 3, "batch_4")),
    split_gene_chromosome = dplyr::mutate(index, chromosome = replace(chromosome, 3, "chr1")),
    other_chromosome = dplyr::mutate(index, chromosome = replace(chromosome, 1, "chr23")),
    chromosome_path_text = dplyr::mutate(index, chromosome = replace(chromosome, 1, "../outside/path")),
    failed_outcome = dplyr::mutate(index, status = replace(status, 1, "FAILED")),
    unknown_status = dplyr::mutate(index, status = replace(status, 1, "SKIPPED")),
    missing_status = dplyr::mutate(index, status = replace(status, 1, NA_character_)),
    inconsistent_convergence = dplyr::mutate(index, converged = replace(converged, 1, FALSE)),
    missing_convergence = dplyr::mutate(index, converged = replace(converged, 1, NA)),
    all_identity_columns_missing = tibble::tibble(other_metadata = c("first", "second")))
  expected_counts <- c(0L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L,
                       3L, 3L, 3L, 3L, 3L, 1L)
  for (i in seq_along(accepted_indices)) {
    candidate <- arguments
    candidate$output_dir <- file.path(root, paste0("accepted ", names(accepted_indices)[[i]]))
    value <- accepted_indices[[i]]
    readr::write_tsv(value, index_file, quote = "needed")
    written <- do.call(writeChromosomeGeneManifests, candidate)
    if (!nrow(value)) {
      stopifnot(is.character(written), length(written) == 0L)
      next
    }
    actual <- bind_rows(lapply(written, read_manifest))
    stopifnot(nrow(actual) == expected_counts[[i]], all(file.exists(written)),
      all(dirname(normalizePath(written)) == normalizePath(candidate$output_dir)),
      all(!grepl("[/\\\\]", basename(written))))
    for (column in c("gene_id", "gene_key", "batch_id")) {
      if (column %in% names(value)) stopifnot(setequal(actual[[column]], value[[column]]))
      else stopifnot(all(is.na(actual[[column]])))
    }
    if (!"gene_key" %in% names(value)) {
      stopifnot(all(is.na(unlist(actual[names(suffixes)], use.names = FALSE))))
    }
  }
  readr::write_tsv(index, index_file)

  # Missing, repeated, foreign, and unmatched artifact locations remain metadata.
  artifact_cases <- list()
  changed <- arguments
  changed$gene_rds_uris <- changed$gene_rds_uris[-1]
  artifact_cases$missing_rds <- changed
  changed <- arguments
  changed$gene_parquet_uris <- changed$gene_parquet_uris[-1]
  artifact_cases$missing_parquet_role <- changed
  changed <- arguments
  changed$gene_rds_uris <- c(changed$gene_rds_uris, changed$gene_rds_uris[[1]])
  artifact_cases$duplicate_uri <- changed
  changed <- arguments
  changed$gene_rds_uris <- c(changed$gene_rds_uris,
    sub("missing-bucket", "another-bucket", changed$gene_rds_uris[[1]]))
  artifact_cases$ambiguous_role <- changed
  changed <- arguments
  changed$gene_parquet_uris[[1]] <- "gs://missing-bucket/gene-foreign.parquet"
  artifact_cases$foreign_gene <- changed
  changed <- arguments
  changed$gene_tsv_uris[[1]] <- "gs://missing-bucket/unrelated.tsv"
  artifact_cases$unknown_suffix <- changed
  changed <- arguments
  changed$gene_log_uris[[1]] <- changed$gene_provenance_uris[[1]]
  changed$gene_provenance_uris[[1]] <- arguments$gene_log_uris[[1]]
  artifact_cases$other_artifact_group <- changed
  for (value in list(NA_character_, "", "line\nlocation", "tab\tlocation", 123, NULL, character())) {
    changed <- arguments
    changed["gene_rds_uris"] <- list(value)
    artifact_cases[[paste0("unmatched_rds_", length(artifact_cases) + 1L)]] <- changed
  }
  for (name in names(artifact_cases)) {
    candidate <- artifact_cases[[name]]
    candidate$output_dir <- file.path(root, paste0("artifact case ", name))
    written <- do.call(writeChromosomeGeneManifests, candidate)
    actual <- bind_rows(lapply(written, read_manifest))
    stopifnot(nrow(actual) == 3L, setequal(actual$gene_id, genes$gene_id))
    if (name %in% c("missing_rds", "ambiguous_role")) {
      stopifnot(is.na(actual$rds_uri[actual$gene_key == genes$gene_key[[3L]]]))
    }
    if (startsWith(name, "unmatched_rds_")) stopifnot(all(is.na(actual$rds_uri)))
  }
  message("Chromosome gene manifest metadata merge checks passed.")
}
