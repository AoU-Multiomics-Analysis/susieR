# Genotype locations are manifest data until WDL constructs typed File inputs.
run_pgen_manifest_tests <- function(repo_dir = ".") {
  suppressPackageStartupMessages(library(dplyr))
  functions_dir <- Sys.getenv("SUSIER_FUNCTIONS_PATH",
                              unset = file.path(repo_dir, "R", "utils"))
  scripts_dir <- Sys.getenv("SUSIER_SCRIPTS_PATH",
                            unset = file.path(repo_dir, "R", "scripts"))
  source(file.path(functions_dir, "PgenBatchFunctions.R"))
  stopifnot("The PGEN manifest reader must exist" = exists("readPgenManifest", mode = "function"),
            "The PGEN manifest writer must exist" = exists("writePgenManifest", mode = "function"))
  root <- tempfile("pgen manifest 'quoted' $literal ")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  manifest_file <- file.path(root, "manifest.tsv")
  normalized_file <- file.path(root, "normalized.tsv")
  manifest <- tibble::tibble(
    chromosome = c("chr2", "1"),
    pgen = c("gs://missing-bucket/chr2 'quoted' $(touch PWNED) `touch PWNED`.pgen",
             "/missing local directory/chr1.pgen"),
    pvar = c('gs://missing-bucket/chr2 "literal quote".pvar', "/missing local directory/chr1.pvar"),
    psam = rep("gs://missing-bucket/shared/NA", 2L))
  emit <- function(value) readr::write_tsv(value, manifest_file, quote = "all")
  read_manifest <- function() readPgenManifest(manifest_file)
  emit(manifest)
  parsed <- read_manifest()
  stopifnot(identical(lapply(parsed, identity), lapply(manifest, identity)))
  writePgenManifest(parsed, normalized_file)
  expected_lines <- apply(as.data.frame(manifest), 1L, paste, collapse = "\t")
  stopifnot(identical(readLines(normalized_file, warn = FALSE), unname(expected_lines)))
  # Only the manifest File may reach the path helper. Genotype locations do not exist.
  guarded_read <- readPgenManifest
  guarded_env <- new.env(parent = environment(guarded_read))
  local_file <- pgenLocalFile
  guarded_env$pgenLocalFile <- function(path, ...) {
    stopifnot("The planner must not open manifest genotype paths" = identical(path, manifest_file))
    local_file(path, ...)
  }
  environment(guarded_read) <- guarded_env
  stopifnot(identical(lapply(guarded_read(manifest_file), identity), lapply(manifest, identity)))
  error <- function(code, pattern) {
    result <- tryCatch({force(code); NULL}, error = identity)
    stopifnot(inherits(result, "error"), grepl(pattern, conditionMessage(result), ignore.case = TRUE))
    invisible(result)
  }
  invalid <- list(
    empty = manifest[FALSE, ],
    missing_column = dplyr::select(manifest, -psam),
    extra_column = dplyr::mutate(manifest, extra = "bad"),
    reordered_columns = dplyr::select(manifest, pgen, chromosome, pvar, psam),
    duplicate_chromosome = dplyr::mutate(manifest, chromosome = c("1", "1")),
    duplicate_alias = dplyr::mutate(manifest, chromosome = c("chr1", "1")),
    unsupported_chromosome = dplyr::mutate(manifest, chromosome = c("X", "1")),
    zero_chromosome = dplyr::mutate(manifest, chromosome = c("0", "1")),
    missing_chromosome = dplyr::mutate(manifest, chromosome = c(NA_character_, "1")),
    missing_path = dplyr::mutate(manifest, pgen = c(NA_character_, pgen[[2]])),
    blank_path = dplyr::mutate(manifest, pgen = c("", pgen[[2]])),
    relative_path = dplyr::mutate(manifest, pgen = c("chr2.pgen", pgen[[2]])),
    tilde_path = dplyr::mutate(manifest, pgen = c("~/chr2.pgen", pgen[[2]])),
    literal_na_path = dplyr::mutate(manifest, pgen = c("NA", pgen[[2]])),
    empty_uri = dplyr::mutate(manifest, pgen = c("gs://", pgen[[2]])),
    tab_path = dplyr::mutate(manifest, pgen = c("gs://bucket/bad\tpath", pgen[[2]])),
    newline_path = dplyr::mutate(manifest, pgen = c("gs://bucket/bad\npath", pgen[[2]])),
    carriage_path = dplyr::mutate(manifest, pgen = c("gs://bucket/bad\rpath", pgen[[2]])))
  for (name in names(invalid)) {
    emit(invalid[[name]])
    error(read_manifest(), "manifest|schema|column|chromosome|autosome|duplicate|nonempty|absolute|single")
    error(writePgenManifest(invalid[[name]], normalized_file),
          "manifest|schema|column|chromosome|autosome|duplicate|nonempty|absolute|single")
  }
  emit(manifest)
  error(writePgenManifest(manifest, "gs://bucket/unlocalized-output"), "local")
  error(readPgenManifest("gs://bucket/unlocalized-manifest"), "localization")
  # A literal NA is a relative location, not a missing field after parsing.
  emit(dplyr::mutate(manifest, pgen = c("NA", pgen[[2]])))
  error(read_manifest(), "absolute")
  # Unequal field counts must not be silently padded or combined.
  writeLines(c("chromosome\tpgen\tpvar\tpsam", "1\tgs://bucket/a\tgs://bucket/b"), manifest_file)
  error(read_manifest(), "manifest|parse|schema|nonempty")
  writeLines(c("chromosome\tpgen\tpvar\tpsam", "1\tgs://bucket/a\tgs://bucket/b\tgs://bucket/c\textra"), manifest_file)
  error(read_manifest(), "manifest|parse|schema|single")
  emit(manifest)

  # Compare the real planner CLI with its existing chromosome-list interface.
  bed <- file.path(root, "phenotypes.bed")
  qtl <- file.path(root, "qtl.tsv")
  chromosomes_file <- file.path(root, "chromosomes.txt")
  writeLines(c("#chr\tstart\tend\tphenotype_id\tS1", "chr1\t100\t101\ta\t1",
               "chr1\t200\t201\tb\t2", "chr2\t300\t301\tc\t3"), bed)
  readr::write_tsv(tibble::tibble(phenotype_id = c("a", "b", "c"), qval = .01), qtl)
  writeLines(manifest$chromosome, chromosomes_file)
  previous <- Sys.getenv("SUSIER_FUNCTIONS_PATH", unset = NA)
  Sys.setenv(SUSIER_FUNCTIONS_PATH = functions_dir)
  on.exit(if (is.na(previous)) Sys.unsetenv("SUSIER_FUNCTIONS_PATH") else
    Sys.setenv(SUSIER_FUNCTIONS_PATH = previous), add = TRUE)
  invoke <- function(label, args, expected_error = NULL) {
    log <- file.path(root, paste0(label, ".log"))
    code <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
      shQuote(c(file.path(scripts_dir, "plan_pgen_susie.R"), "--phenotype-bed", bed,
        "--qtl-results", qtl, "--cis-distance", "10", args)), stdout = log, stderr = log))
    if (is.null(expected_error) && !identical(code, 0L)) {
      stop("Planner CLI failed: ", paste(tail(readLines(log, warn = FALSE), 20L), collapse = "\n"), call. = FALSE)
    }
    if (!is.null(expected_error)) {
      stopifnot(code != 0L, grepl(expected_error, paste(readLines(log, warn = FALSE), collapse = "\n"),
        ignore.case = TRUE))
    }
  }
  plan_manifest <- file.path(root, "manifest plan")
  plan_old <- file.path(root, "old plan")
  invoke("manifest", c("--pgen-manifest", manifest_file, "--output-dir", plan_manifest))
  invoke("old", c("--chromosomes-file", chromosomes_file, "--output-dir", plan_old))
  stopifnot(identical(readLines(file.path(plan_manifest, "chromosome_pgens.tsv")), unname(expected_lines)))
  for (name in c("expected_outcomes.tsv", "batch_chromosome_indices.json", "batch_chromosomes.txt")) {
    stopifnot(identical(readLines(file.path(plan_manifest, name)), readLines(file.path(plan_old, name))))
  }
  indices <- jsonlite::read_json(file.path(plan_manifest, "batch_chromosome_indices.json"), simplifyVector = TRUE)
  stopifnot(identical(indices, c(1L, 0L)))
  for (i in seq_len(2L)) {
    name <- sprintf("batch_%06d", i)
    stopifnot(identical(readLines(file.path(plan_manifest, "batches", paste0(name, ".tsv"))),
      readLines(file.path(plan_old, "batches", paste0(name, ".tsv")))),
      identical(readLines(file.path(plan_manifest, "batches", paste0(name, ".bed"))),
      readLines(file.path(plan_old, "batches", paste0(name, ".bed")))))
  }
  invoke("neither", c("--output-dir", file.path(root, "neither")), "exactly one")
  invoke("both", c("--pgen-manifest", manifest_file, "--chromosomes-file", chromosomes_file,
    "--output-dir", file.path(root, "both")), "exactly one")
  invoke("unlocalized", c("--pgen-manifest", "gs://bucket/unlocalized-manifest",
    "--output-dir", file.path(root, "unlocalized")), "localization")
  stopifnot(!file.exists("PWNED"))
  message("PGEN manifest parser and planner CLI tests passed: typed-path data, literal URIs, and original plan parity")
  invisible(TRUE)
}
