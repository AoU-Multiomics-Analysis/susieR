# Check batch planning and preserve input rows during aggregation.
run_pgen_batch_tests <- function(repo_dir) {
  suppressPackageStartupMessages(library(dplyr))
  functions_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH",
                               unset = file.path(repo_dir, "R", "utils"))
  scripts_path <- Sys.getenv("SUSIER_SCRIPTS_PATH",
                             unset = file.path(repo_dir, "R", "scripts"))
  source(file.path(functions_path, "ImportFunctions.R"))
  source(file.path(functions_path, "InitFunctions.R"))
  source(file.path(functions_path, "PgenBatchFunctions.R"))
  stopifnot(exists("planPgenBatches"), exists("writePgenBatchRegions"),
            exists("aggregatePgenBatches"))
  root <- tempfile("pgen batches ")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  old <- setwd(root)
  on.exit(setwd(old), add = TRUE)
  error <- function(expr, pattern = NULL) {
    result <- tryCatch({
      force(expr)
      NULL
    }, error = identity)
    stopifnot(inherits(result, "error"))
    if (!is.null(pattern)) stopifnot(grepl(pattern, conditionMessage(result)))
    invisible(result)
  }
  read_tsv <- function(path) readr::read_tsv(path, show_col_types = FALSE)
  ids <- sprintf("p%03d", 101:1)
  bed_lines <- c("#chr\tstart\tend\tphenotype_id\tS 2\tS1",
                 paste("chr1", seq_along(ids) * 100, seq_along(ids) * 100 + 1,
                       ids, "1.00000000", "NA", sep = "\t"))
  writeLines(bed_lines, "phenotypes file.bed")
  readr::write_tsv(tibble(phenotype_id = ids, qval = .001), "qtl file.tsv")
  plan <- planPgenBatches("phenotypes file.bed", "qtl file.tsv", "chr1", 200,
                         "plan")
  expected <- read_tsv("plan/expected_outcomes.tsv")
  stopifnot(identical(as.integer(table(expected$batch_id)), c(50L, 50L, 1L)),
            identical(expected$phenotype_id, ids),
            expected$phenotype_pos[[1]] == 100,
            expected$window_start[[1]] == 1,
            expected$requested_window_start[[1]] == -100,
            expected$region[[1]] == "chr1:-100-300")
  paths <- readLines("plan/batch_phenotype_beds.txt")
  worklists <- readLines("plan/batch_worklists.txt")
  stopifnot(length(paths) == 3L, all(file.exists(paths)), all(file.exists(worklists)),
            identical(readLines(paths[[1]]), bed_lines[1:51]),
            identical(readLines(paths[[3]]), bed_lines[c(1, 102)]),
            identical(readLines("plan/batch_chromosomes.txt"), rep("chr1", 3)))
  planPgenBatches("phenotypes file.bed", "qtl file.tsv", "chr1", 200,
                  "plan100", batch_size = 100)
  stopifnot(identical(as.integer(table(read_tsv("plan100/expected_outcomes.tsv")$batch_id)),
                      c(100L, 1L)))
  small <- c("#chr\tstart\tend\tID\tS1", "chr2\t100\t100\tz\t1",
             "chr1\t20\t22\tb\t2", "chr1\t10\t12\ta\t3",
             "chr1\t1000\t1000\tc\t4", "chr1\t20\t22\td\t5")
  writeLines(small, "small.bed")
  readr::write_tsv(tibble(phenotype_id = c("z", "b", "a", "c", "d"),
                         qval = c(.01, .01, .01, .01, .9)), "small.tsv")
  readr::write_tsv(tibble(phenotype_id = c("z", "b", "a", "c", "d"),
                         gene_id = c("Z", "A", "A", "C", "D")), "map.tsv")
  planPgenBatches("small.bed", "small.tsv", c("chr2", "chr1"), 5,
                  "small plan", batch_size = 1, gene_map = "map.tsv")
  small_expected <- read_tsv("small plan/expected_outcomes.tsv")
  stopifnot(identical(small_expected$phenotype_id, c("a", "b", "c", "z")),
            identical(small_expected$batch_id, c("batch_000001", "batch_000001",
                                               "batch_000002", "batch_000003")),
            small_expected$phenotype_pos[[1]] == 11)
  writePgenBatchRegions(small_expected, "regions.bed1")
  regions <- readr::read_tsv("regions.bed1", col_names = FALSE, show_col_types = FALSE)
  stopifnot(nrow(regions) == 4, identical(regions$X2, c(6, 16, 995, 95)),
            identical(regions$X3, c(16, 26, 1005, 105)),
            !anyDuplicated(regions$X4), all(grepl("^[A-Za-z0-9_]+$", regions$X4)))
  error(planPgenBatches("small.bed", "small.tsv", "chr1", 5, "bad"), "chromosome")
  readr::write_tsv(tibble(phenotype_id = "a", gene_id = "A"), "missing map.tsv")
  error(planPgenBatches("small.bed", "small.tsv", c("chr1", "chr2"), 5,
                       "bad", gene_map = "missing map.tsv"), "mapping")
  readr::write_tsv(tibble(phenotype_id = c("z", "b", "a", "c", "d"),
                         gene_id = c("A", "A", "A", "C", "D")), "cross map.tsv")
  error(planPgenBatches("small.bed", "small.tsv", c("chr1", "chr2"), 5,
                       "bad", gene_map = "cross map.tsv"), "chromosome")
  writeLines(c(small, small[[3]]), "duplicate.bed")
  error(planPgenBatches("duplicate.bed", "small.tsv", c("chr1", "chr2"), 5,
                       "bad"), "duplicate")
  cluster_ids <- c("x:clu_1:a", "x:clu_1:b", "x:clu_2:c")
  writeLines(c(small[[1]], paste("chr1", 10:12, 11:13, cluster_ids, 1, sep = "\t")),
             "cluster.bed")
  readr::write_tsv(tibble(phenotype_id = cluster_ids, qval = c(.01, .02, .01),
                         pval_beta = c(.02, .001, .03)), "cluster.tsv")
  planPgenBatches("cluster.bed", "cluster.tsv", "chr1", 5, "cluster plan",
                  select_top_phenotype_per_cluster = TRUE,
                  top_phenotype_pvalue_column = "pval_beta")
  stopifnot(identical(read_tsv("cluster plan/expected_outcomes.tsv")$phenotype_id,
                      cluster_ids[c(2, 3)]))

  # Numeric chromosome order and same-position ties must not depend on BED order.
  writeLines(c(small[[1]], "chr10\t10\t12\tp10\t1", "chr2\t10\t12\tp2z\t2",
               "chr2\t10\t12\tp2a\t3"), "ties.bed")
  readr::write_tsv(tibble(phenotype_id = c("p10", "p2z", "p2a"), qval = .01), "ties.tsv")
  planPgenBatches("ties.bed", "ties.tsv", c("chr10", "chr2"), 5, "ties plan", batch_size = 50)
  ties <- read_tsv("ties plan/expected_outcomes.tsv")
  stopifnot(identical(ties$phenotype_id, c("p2a", "p2z", "p10")),
            identical(ties$batch_id, c("batch_000001", "batch_000001", "batch_000002")))
  overlap <- tibble(chromosome = "chr1", window_start = c(1, 5, 100), window_end = c(10, 20, 110))
  writePgenBatchRegions(overlap, "overlap.bed1")
  stopifnot(identical(readLines("overlap.bed1"),
                     c("chr1\t1\t10\tregion_000000001", "chr1\t5\t20\tregion_000000002",
                       "chr1\t100\t110\tregion_000000003")))
  error(planPgenBatches("small.bed", "small.tsv", c("chr1", "chr2"), 5,
                       "bad", batch_size = 0), "batch_size")
  error(planPgenBatches("small.bed", "small.tsv", c("chr1", "chr2"), -1,
                       "bad"), "cis_distance")
  error(writePgenBatchRegions(mutate(overlap, window_start = 0), "bad.bed1"), "window")
  # Chunk splitting must preserve rows on both sides of the chunk boundary.
  many_ids <- paste0("m", seq_len(1201))
  many_lines <- c(small[[1]], paste("chr1", 1:1201, 2:1202, many_ids, "1e-10", sep = "\t"))
  writeLines(many_lines, "many.bed")
  readr::write_tsv(tibble(phenotype_id = many_ids[c(1000, 1001, 1201)], qval = .01), "many.tsv")
  planPgenBatches("many.bed", "many.tsv", "chr1", 5, "many plan")
  stopifnot(identical(readLines(readLines("many plan/batch_phenotype_beds.txt")),
                     many_lines[c(1, 1001, 1002, 1202)]))

  compressed <- gzfile("small.bed.gz", "wt")
  writeLines(small, compressed)
  close(compressed)
  planPgenBatches("small.bed.gz", "small.tsv", c("chr1", "chr2"), 5, "gz plan")
  gz_output <- readLines("gz plan/batch_phenotype_beds.txt")
  stopifnot(identical(readLines(gz_output[[1]]), small[c(1, 3, 4, 5)]))

  # Standard compressed phenotype BEDs also use the .bed.bgz suffix.
  Rsamtools::bgzip("small.bed", dest = "small.bed.bgz", overwrite = TRUE)
  planPgenBatches("small.bed.bgz", "small.tsv", c("chr1", "chr2"), 5, "bgz plan")
  bgz_output <- readLines("bgz plan/batch_phenotype_beds.txt")
  stopifnot(identical(readLines(bgz_output[[1]]), small[c(1, 3, 4, 5)]),
            identical(read_tsv("bgz plan/expected_outcomes.tsv")$phenotype_pos,
                      c(11, 21, 1000, 100)))
  # Use exact chromosome map keys while keeping legacy BED chromosome in regions.
  planPgenBatches("phenotypes file.bed", "qtl file.tsv", "1", 200, "numeric key plan")
  numeric_key <- read_tsv("numeric key plan/expected_outcomes.tsv")
  stopifnot(all(numeric_key$chromosome == "1"),
            numeric_key$region[[1]] == "chr1:-100-300",
            identical(readLines("numeric key plan/batch_chromosomes.txt"), rep("1", 3)))
  writeLines(sub("^chr1\t", "1\t", bed_lines), "numeric chromosome.bed")
  planPgenBatches("numeric chromosome.bed", "qtl file.tsv", "chr1", 200, "prefixed key plan")
  prefixed_key <- read_tsv("prefixed key plan/expected_outcomes.tsv")
  stopifnot(all(prefixed_key$chromosome == "chr1"), prefixed_key$region[[1]] == "1:-100-300")

  # The expected-outcomes path is a compatibility input, not a validation plan.
  writeLines("Readable compatibility input", "expected.tsv")
  wanted <- small_expected %>% filter(batch_id %in% c("batch_000001", "batch_000002"))
  index <- wanted %>% transmute(batch_id, chromosome, gene_id, gene_key = paste0("gene-", gene_id),
    outcome_id, phenotype_id, seed = 1, status = c("CONVERGED", "NONCONVERGED", "CONVERGED"),
    converged = c(TRUE, FALSE, TRUE), credible_set_count = c(1, 0, 0), reason = "")
  own <- wanted %>% select(batch_id, gene_id, outcome_id, phenotype_id)
  summary <- bind_cols(own[1, ], tibble(molecular_trait_id = "a", variant = "1_10_A_G",
    chromosome = "1", position = 10L, ref = "A", alt = "G", cs_id = "a_1", cs_index = "1",
    region = "chr1:6-16", pip = .9, z = 4, cs_min_r2 = .8, cs_avg_r2 = .9, cs_size = 1L,
    posterior_mean = .4, posterior_sd = .1, cs_log10bf = 3,
    has_covariate_explained_variant = FALSE, phenotype_residual_SD = 1.25))
  cs <- bind_cols(own[1, ], tibble(molecular_trait_id = "a", cs_id = "a_1", cs_index = "1",
    region = "chr1:6-16", cs_log10bf = 3, cs_avg_r2 = .9, cs_min_r2 = .8, cs_size = 1,
    low_purity = FALSE, has_covariate_explained_variant = FALSE))
  diagnostic <- own[1, ] %>% mutate(molecular_trait_id = "a", variant = "1_11_A_T",
    input_scale = 1, residual_sd = 0, residual_sd_ratio = 0)
  summary_tables <- list(
    bind_rows(summary, summary, mutate(summary, pip = .2, z = -2, posterior_mean = -.3,
                              phenotype_residual_SD = 2.5)) %>%
      mutate(feature_converged = c(TRUE, TRUE, FALSE), first_note = c("first", "first", "different")),
    mutate(summary, molecular_trait_id = "foreign", outcome_id = "outside",
      phenotype_id = "foreign", gene_id = "unlisted gene", cs_id = "unknown set",
      cs_size = 99L, feature_converged = TRUE, second_note = "second"))
  cs_tables <- list(
    bind_rows(cs, cs, mutate(cs, cs_size = -5, low_purity = TRUE)) %>%
      mutate(feature_converged = c(FALSE, FALSE, TRUE)),
    mutate(cs, molecular_trait_id = "foreign", outcome_id = "outside",
      cs_id = "unknown set", feature_converged = TRUE))
  index_tables <- list(
    index[c(2, 1), ] %>% mutate(status = "FAILED", converged = TRUE,
                               credible_set_count = 77),
    bind_rows(index[3, ], index[3, ]))
  diagnostic_tables <- list(bind_rows(diagnostic, diagnostic,
      mutate(diagnostic, input_scale = 2, residual_sd = .1, residual_sd_ratio = .05)),
    mutate(diagnostic, gene_id = "foreign", input_scale = -1))
  emit <- function(sums = summary_tables, sets = cs_tables, indices = index_tables,
                   diagnostics = diagnostic_tables) {
    for (i in 1:2) {
      arrow::write_parquet(sums[[i]], paste0(i, " summary.parquet"))
      arrow::write_parquet(sets[[i]], paste0(i, " cs.parquet"))
      readr::write_tsv(indices[[i]], paste0(i, " index.tsv"), quote = "needed")
      readr::write_tsv(diagnostics[[i]], paste0(i, " diagnostic.tsv"), quote = "needed")
    }
  }
  aggregate <- function(summary_paths = c("1 summary.parquet", "2 summary.parquet"),
                        cs_paths = c("1 cs.parquet", "2 cs.parquet"),
                        index_paths = c("1 index.tsv", "2 index.tsv"),
                        diagnostic_paths = c("1 diagnostic.tsv", "2 diagnostic.tsv"),
                        expected_path = "expected.tsv") {
    aggregatePgenBatches(expected_path, summary_paths, cs_paths, index_paths,
                         diagnostic_paths, "aggregate")
  }
  without_warnings <- function(expr) {
    warnings <- character()
    value <- withCallingHandlers(force(expr), warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    })
    stopifnot("Loading and merging must not warn about scientific metadata" = !length(warnings))
    value
  }
  check_tsv <- function(path, expected) {
    column_types <- do.call(readr::cols, lapply(expected, function(value) {
      switch(typeof(value), double = readr::col_double(), integer = readr::col_integer(),
        logical = readr::col_logical(), readr::col_character())
    }))
    actual <- tibble::as_tibble(readr::read_tsv(path, col_types = column_types,
      na = "NA", show_col_types = FALSE))
    stopifnot(identical(as.list(actual), as.list(expected)))
  }
  emit()
  without_warnings(aggregate())
  # Repeated and foreign rows, differing values, and inconsistent metadata stay intact.
  stopifnot(identical(arrow::read_parquet("aggregate/susie.parquet"), bind_rows(summary_tables)),
    identical(arrow::read_parquet("aggregate/credible_sets.parquet"), bind_rows(cs_tables)))
  check_tsv("aggregate/outcome_index.tsv", bind_rows(index_tables))
  check_tsv("aggregate/covariate_explained_variants.tsv", bind_rows(diagnostic_tables))

  # Each list controls its own row order. A repeated path contributes its rows again.
  without_warnings(aggregate(summary_paths = c("2 summary.parquet", "1 summary.parquet", "1 summary.parquet"),
    cs_paths = c("2 cs.parquet", "1 cs.parquet", "1 cs.parquet"),
    index_paths = c("2 index.tsv", "1 index.tsv", "1 index.tsv"),
    diagnostic_paths = c("2 diagnostic.tsv", "1 diagnostic.tsv", "1 diagnostic.tsv")))
  stopifnot(identical(arrow::read_parquet("aggregate/susie.parquet"),
                      bind_rows(summary_tables[c(2, 1, 1)])),
    identical(arrow::read_parquet("aggregate/credible_sets.parquet"),
              bind_rows(cs_tables[c(2, 1, 1)])))
  check_tsv("aggregate/outcome_index.tsv", bind_rows(index_tables[c(2, 1, 1)]))
  check_tsv("aggregate/covariate_explained_variants.tsv", bind_rows(diagnostic_tables[c(2, 1, 1)]))
  # Expected-outcome content does not set row order or enforce completeness.
  writeLines(character(), "expected.tsv")
  without_warnings(aggregate())
  check_tsv("aggregate/outcome_index.tsv", bind_rows(index_tables))

  # Older files receive unknown SD values. Supplied SD values keep their type.
  mixed_sd <- list(select(summary_tables[[1L]], -phenotype_residual_SD), summary_tables[[2L]])
  emit(sums = mixed_sd)
  without_warnings(aggregate())
  stopifnot(identical(arrow::read_parquet("aggregate/susie.parquet"), bind_rows(mixed_sd)))
  text_sd <- lapply(summary_tables, function(table) {
    mutate(table, phenotype_residual_SD = as.character(phenotype_residual_SD))
  })
  emit(sums = text_sd)
  without_warnings(aggregate())
  stopifnot(identical(arrow::read_parquet("aggregate/susie.parquet"), bind_rows(text_sd)))
  mixed_text_sd <- list(select(summary_tables[[1L]], -phenotype_residual_SD), text_sd[[2L]])
  emit(sums = mixed_text_sd)
  without_warnings(aggregate())
  stopifnot(identical(arrow::read_parquet("aggregate/susie.parquet"), bind_rows(mixed_text_sd)))

  # Fill absent metadata only. Conflicting or unknown index values stay unknown.
  metadata_index <- tibble(outcome_id = c("a", "a", "b", "conflict", "conflict"),
    converged = c(TRUE, TRUE, FALSE, TRUE, FALSE),
    status = c("FAILED", "FAILED", "CONVERGED", "FAILED", "NONCONVERGED"))
  unknown_index <- tibble(outcome_id = "without convergence", status = "FAILED")
  metadata_summary <- tibble(molecular_trait_id = c("a", "a", "b", "conflict", "unknown",
                                                    "without convergence"),
    variant = "same label", pip = c(.1, .2, .3, .4, .5, .6),
    feature_converged = c(NA, FALSE, NA, NA, NA, NA))
  metadata_summary2 <- tibble(molecular_trait_id = "b", variant = "other label", pip = .7)
  metadata_cs <- tibble(molecular_trait_id = c("a", "b", "conflict", "unknown", "without convergence"),
    cs_id = c("a", "b", "c", "d", "e"))
  metadata_cs2 <- tibble(molecular_trait_id = "a", cs_id = "explicit", feature_converged = FALSE)
  emit(sums = list(metadata_summary, metadata_summary2), sets = list(metadata_cs, metadata_cs2),
    indices = list(metadata_index, unknown_index))
  without_warnings(aggregate())
  actual <- arrow::read_parquet("aggregate/susie.parquet")
  expected <- bind_rows(metadata_summary, metadata_summary2) %>%
    mutate(feature_converged = c(TRUE, FALSE, FALSE, NA, NA, NA, FALSE),
           phenotype_residual_SD = NA_real_)
  stopifnot(identical(actual[names(expected)], expected),
    is.double(actual$phenotype_residual_SD), all(is.na(actual$phenotype_residual_SD)))
  actual_cs <- arrow::read_parquet("aggregate/credible_sets.parquet")
  expected_cs <- bind_rows(metadata_cs, metadata_cs2) %>%
    mutate(feature_converged = c(TRUE, FALSE, NA, NA, NA, FALSE))
  stopifnot(identical(actual_cs[names(expected_cs)], expected_cs))
  check_tsv("aggregate/outcome_index.tsv", bind_rows(metadata_index, unknown_index))

  # Empty tables do not require successful statuses or credible-set counts.
  empty_summary <- tibble(molecular_trait_id = character(), variant = character())
  empty_cs <- tibble(molecular_trait_id = character(), cs_id = character())
  emit(sums = list(empty_summary, empty_summary), sets = list(empty_cs, empty_cs))
  without_warnings(aggregate())
  empty <- arrow::read_parquet("aggregate/susie.parquet")
  stopifnot(nrow(empty) == 0L, is.double(empty$phenotype_residual_SD),
    is.logical(empty$feature_converged),
    nrow(arrow::read_parquet("aggregate/credible_sets.parquet")) == 0L)
  check_tsv("aggregate/outcome_index.tsv", bind_rows(index_tables))

  # File readability and WDL localization errors still stop loading.
  emit()
  error(aggregate(expected_path = "missing expected.tsv"), "read|exist|local")
  error(aggregate(expected_path = "gs://bucket/expected.tsv"), "localization")
  for (argument in c("summary_paths", "cs_paths", "index_paths", "diagnostic_paths")) {
    invalid <- setNames(list("missing batch input"), argument)
    error(do.call(aggregate, invalid), "read|exist|local")
    invalid[[1L]] <- "gs://bucket/batch-input"
    error(do.call(aggregate, invalid), "localization")
  }
  error(planPgenBatches("gs://bucket/phenotype.bed", "small.tsv", "chr1", 5,
                       "bad"), "localization")

  # The command entrypoints preserve named paths, including spaces and repeats.
  chr_path <- file.path(root, "chromosomes file.txt")
  writeLines("chr1", chr_path)
  previous_functions <- Sys.getenv("SUSIER_FUNCTIONS_PATH", unset = NA_character_)
  Sys.setenv(SUSIER_FUNCTIONS_PATH = functions_path)
  on.exit(if (is.na(previous_functions)) Sys.unsetenv("SUSIER_FUNCTIONS_PATH") else
    Sys.setenv(SUSIER_FUNCTIONS_PATH = previous_functions), add = TRUE)
  run_cli <- function(script, args) {
    output <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                        c(shQuote(file.path(scripts_path, script)), shQuote(args)),
                        stdout = TRUE, stderr = TRUE))
    status <- attr(output, "status")
    if (!is.null(status) && status != 0L) stop(paste(output, collapse = "\n"))
    output
  }
  run_cli("plan_pgen_susie.R", c("--phenotype-bed", "phenotypes file.bed", "--qtl-results",
                                "qtl file.tsv", "--chromosomes-file", chr_path, "--cis-distance", "200",
                                "--output-dir", "cli plan"))
  stopifnot(nrow(read_tsv("cli plan/expected_outcomes.tsv")) == 101L)
  emit()
  writeLines(c("2 summary.parquet", "1 summary.parquet", "1 summary.parquet"), "summaries list.txt")
  writeLines(c("2 cs.parquet", "1 cs.parquet", "1 cs.parquet"), "sets list.txt")
  writeLines(c("2 index.tsv", "1 index.tsv", "1 index.tsv"), "indices list.txt")
  writeLines(c("2 diagnostic.tsv", "1 diagnostic.tsv", "1 diagnostic.tsv"), "diagnostics list.txt")
  run_cli("aggregate_pgen_susie.R", c("--expected-outcomes", "expected.tsv",
      "--summary-files", "summaries list.txt", "--credible-set-files", "sets list.txt",
      "--index-files", "indices list.txt", "--diagnostic-files", "diagnostics list.txt",
      "--output-dir", "cli aggregate"))
  stopifnot(identical(arrow::read_parquet("cli aggregate/susie.parquet"),
                      bind_rows(summary_tables[c(2, 1, 1)])))
  check_tsv("cli aggregate/outcome_index.tsv", bind_rows(index_tables[c(2, 1, 1)]))
  message("PGEN batch planning and simple aggregation tests passed")
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  run_pgen_batch_tests(normalizePath(if (length(args)) args[[1L]] else ".", mustWork = TRUE))
}
