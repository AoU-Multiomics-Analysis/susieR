# These fixtures catch split genes, lost cis windows and incomplete aggregation.
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

  # Each file set has one file per batch, paired in the supplied vector order.
  wanted <- small_expected %>% filter(batch_id %in% c("batch_000001", "batch_000002"))
  readr::write_tsv(wanted, "expected.tsv")
  index <- wanted %>% transmute(batch_id, chromosome, gene_id, gene_key = paste0("gene-", gene_id),
                               outcome_id, phenotype_id, seed = 1L,
                               status = c("CONVERGED", "NONCONVERGED", "CONVERGED"),
                               converged = c(TRUE, FALSE, TRUE),
                               credible_set_count = c(1L, 0L, 0L), reason = "")
  own <- wanted %>% select(batch_id, gene_id, outcome_id, phenotype_id)
  summary <- tibble(molecular_trait_id = "a", variant = "1_10_A_G", chromosome = "1",
    position = 10L, ref = "A", alt = "G", cs_id = "a_1", cs_index = "1", region = "chr1:6-16",
    pip = .9, z = 4, cs_min_r2 = .8, cs_avg_r2 = .9, cs_size = 1L,
    posterior_mean = .4, posterior_sd = .1, cs_log10bf = 3,
    has_covariate_explained_variant = FALSE, phenotype_residual_SD = 1.25)
  summary <- bind_cols(own[1, ], summary)
  cs <- tibble(molecular_trait_id = "a", cs_id = "a_1", cs_index = 1,
               region = "chr1:6-16", cs_log10bf = 3, cs_avg_r2 = .9, cs_min_r2 = .8,
               cs_size = 1, low_purity = FALSE, has_covariate_explained_variant = FALSE)
  cs <- bind_cols(own[1, ], cs)
  diagnostic <- own[1, ] %>% mutate(molecular_trait_id = "a", variant = "1_11_A_T",
                                    input_scale = 1, residual_sd = 0,
                                    residual_sd_ratio = 0)
  emit <- function(idx = index, sum = summary, sets = cs, diag = diagnostic) {
    for (i in 1:2) {
      batch <- paste0("batch_00000", i)
      arrow::write_parquet(filter(sum, batch_id == batch), paste0(i, " summary.parquet"))
      arrow::write_parquet(filter(sets, batch_id == batch), paste0(i, " cs.parquet"))
      readr::write_tsv(filter(idx, batch_id == batch), paste0(i, " index.tsv"), quote = "needed")
      readr::write_tsv(filter(diag, batch_id == batch), paste0(i, " diagnostic.tsv"), quote = "needed")
    }
  }
  aggregate <- function(...) aggregatePgenBatches(
    "expected.tsv", c("1 summary.parquet", "2 summary.parquet"),
    c("1 cs.parquet", "2 cs.parquet"), c("1 index.tsv", "2 index.tsv"),
    c("1 diagnostic.tsv", "2 diagnostic.tsv"), "aggregate", ...)
  capture_warnings <- function(expr) {
    warnings <- character()
    withCallingHandlers(force(expr), warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    })
    warnings
  }
  emit()
  aggregate()
  stopifnot(nrow(arrow::read_parquet("aggregate/susie.parquet")) == 1,
            nrow(arrow::read_parquet("aggregate/credible_sets.parquet")) == 1,
            identical(read_tsv("aggregate/outcome_index.tsv")$status, index$status),
            nrow(read_tsv("aggregate/covariate_explained_variants.tsv")) == 1,
            identical(arrow::read_parquet("aggregate/susie.parquet")$phenotype_residual_SD, 1.25))
  # A label is not a physical predictor ID. Keep equal and differing rows in order.
  differing_row <- mutate(summary, pip = .2, z = -2, posterior_mean = -.3,
                          posterior_sd = .2)
  repeated_cases <- list(
    identical_rows = bind_rows(summary, summary),
    differing_values = bind_rows(summary, differing_row),
    two_physical_members = mutate(bind_rows(summary, differing_row), cs_size = 2L),
    legacy_join_expansion = mutate(bind_rows(summary, summary, differing_row,
                                             differing_row), cs_size = 2L)
  )
  for (case in names(repeated_cases)) {
    repeated <- repeated_cases[[case]]
    emit(sum = repeated, sets = mutate(cs, cs_size = repeated$cs_size[[1L]]))
    warnings <- capture_warnings(aggregate())
    retained <- arrow::read_parquet("aggregate/susie.parquet")
    stopifnot(identical(retained[names(repeated)], repeated),
      any(grepl("1 summary.parquet", warnings, fixed = TRUE)),
      any(grepl("batch_000001", warnings, fixed = TRUE)),
      any(grepl("outcome=a", warnings, fixed = TRUE)),
      any(grepl(paste0("repeated_row_count=", nrow(repeated) - 1L), warnings, fixed = TRUE)),
      any(grepl("physical member count", warnings, fixed = TRUE)))
  }
  repeated_diagnostic <- bind_rows(diagnostic, diagnostic,
    mutate(diagnostic, input_scale = 2, residual_sd = .1, residual_sd_ratio = .05))
  emit(diag = repeated_diagnostic)
  warnings <- capture_warnings(aggregate())
  stopifnot(identical(tibble::as_tibble(read_tsv("aggregate/covariate_explained_variants.tsv")),
                      repeated_diagnostic),
    any(grepl("1 diagnostic.tsv", warnings, fixed = TRUE)),
    any(grepl("outcome=a", warnings, fixed = TRUE)),
    any(grepl("repeated_row_count=2", warnings, fixed = TRUE)))
  # Long labels must stay complete in data and remain short in warning examples.
  long_variant <- paste0("1_10_A_G_", strrep("long", 300L))
  long_rows <- mutate(bind_rows(summary, summary), variant = long_variant)
  emit(sum = long_rows)
  warnings <- capture_warnings(aggregate())
  stopifnot(identical(arrow::read_parquet("aggregate/susie.parquet")$variant,
                      rep(long_variant, 2L)),
    !any(grepl(long_variant, warnings, fixed = TRUE)),
    max(nchar(warnings)) < 700L)
  # Repeated labels cannot hide a definite excess or shortage of CS members.
  emit(sum = bind_rows(summary, summary, mutate(summary, variant = "1_12_A_T",
                                                position = 12L)))
  suppressWarnings(error(aggregate(), "missing|excess"))
  emit(sum = mutate(bind_rows(summary, differing_row), cs_size = 3L),
       sets = mutate(cs, cs_size = 3L))
  suppressWarnings(error(aggregate(), "missing|excess"))
  emit(idx = mutate(index, credible_set_count = if_else(outcome_id == "a", 2L, credible_set_count)),
       sets = bind_rows(cs, cs))
  error(aggregate(), "Duplicate credible set")
  emit()
  # Backfill old files from the index, including retained NONCONVERGED sets.
  convergence_index <- index %>%
    mutate(status = if_else(outcome_id == "a", "NONCONVERGED", "CONVERGED"),
           converged = outcome_id != "a", credible_set_count = 1L) %>%
    slice(3:1)
  convergence_summary <- bind_rows(summary,
    bind_cols(own[2L, ], select(summary, -all_of(names(own)))) %>%
      mutate(molecular_trait_id = "b", cs_id = "b_1"),
    bind_cols(own[3L, ], select(summary, -all_of(names(own)))) %>%
      mutate(molecular_trait_id = "c", cs_id = "c_1"))
  convergence_cs <- bind_rows(cs,
    bind_cols(own[2L, ], select(cs, -all_of(names(own)))) %>%
      mutate(molecular_trait_id = "b", cs_id = "b_1"),
    bind_cols(own[3L, ], select(cs, -all_of(names(own)))) %>%
      mutate(molecular_trait_id = "c", cs_id = "c_1"))
  check_convergence <- function() {
    actual_summary <- arrow::read_parquet("aggregate/susie.parquet")
    actual_cs <- arrow::read_parquet("aggregate/credible_sets.parquet")
    stopifnot(identical(actual_summary$feature_converged, c(FALSE, TRUE, TRUE)),
      identical(actual_cs$feature_converged, c(FALSE, TRUE, TRUE)),
      identical(actual_summary[names(convergence_summary)], convergence_summary),
      identical(actual_cs[names(convergence_cs)], convergence_cs))
  }
  emit(idx = convergence_index, sum = convergence_summary, sets = convergence_cs)
  aggregate()
  check_convergence()
  # Old and new files can be combined. Logical NA uses the known index value.
  emit(idx = convergence_index,
    sum = mutate(convergence_summary, feature_converged = c(FALSE, TRUE, NA)),
    sets = mutate(convergence_cs, feature_converged = c(NA, TRUE, TRUE)))
  arrow::write_parquet(filter(convergence_summary, batch_id == "batch_000001"), "1 summary.parquet")
  arrow::write_parquet(filter(convergence_cs, batch_id == "batch_000002"), "2 cs.parquet")
  aggregate()
  check_convergence()
  for (kind in c("summary", "sets")) {
    for (invalid_flag in list(TRUE, "FALSE", 0)) {
      arguments <- list(idx = convergence_index, sum = convergence_summary, sets = convergence_cs)
      arguments[[if (kind == "summary") "sum" else "sets"]] <-
        mutate(arguments[[if (kind == "summary") "sum" else "sets"]], feature_converged = invalid_flag)
      do.call(emit, arguments)
      error(aggregate(), "feature_converged|logical|boolean")
    }
  }
  emit()
  for (kind in c("summary", "cs")) {
    empty <- if (kind == "summary") summary[FALSE, ] else cs[FALSE, ]
    arrow::write_parquet(mutate(empty, feature_converged = character()), paste0("2 ", kind, ".parquet"))
    error(aggregate(), "feature_converged|logical|boolean")
    emit()
  }
  # Missing legacy metadata stays missing, including when all files are old.
  legacy_summary <- select(summary, -phenotype_residual_SD)
  emit(sum = legacy_summary)
  aggregate()
  legacy <- arrow::read_parquet("aggregate/susie.parquet")
  stopifnot("phenotype_residual_SD" %in% names(legacy),
    is.double(legacy$phenotype_residual_SD), all(is.na(legacy$phenotype_residual_SD)))
  emit()
  arrow::write_parquet(filter(legacy_summary, batch_id == "batch_000002"), "2 summary.parquet")
  aggregate()
  stopifnot(identical(arrow::read_parquet("aggregate/susie.parquet")$phenotype_residual_SD, 1.25))
  # Both schemas can have rows; an old outcome must not inherit a new outcome's SD.
  other_summary <- bind_cols(own[3L, ], select(summary, -all_of(names(own)))) %>%
    mutate(molecular_trait_id = "c", variant = "1_20_A_G", position = 20L,
           cs_id = "c_1", phenotype_residual_SD = 2.5)
  other_cs <- bind_cols(own[3L, ], select(cs, -all_of(names(own)))) %>%
    mutate(molecular_trait_id = "c", cs_id = "c_1")
  emit(idx = mutate(index, credible_set_count = if_else(outcome_id == "c", 1L, credible_set_count)),
       sum = bind_rows(summary, other_summary), sets = bind_rows(cs, other_cs))
  arrow::write_parquet(legacy_summary, "1 summary.parquet")
  aggregate()
  stopifnot(identical(arrow::read_parquet("aggregate/susie.parquet")$phenotype_residual_SD,
                      c(NA_real_, 2.5)))
  for (invalid_scale in list("1.25", Inf, -1)) {
    emit(sum = mutate(summary, phenotype_residual_SD = invalid_scale))
    error(aggregate(), "phenotype_residual_SD|numeric|finite|nonnegative")
  }
  emit()
  arrow::write_parquet(mutate(summary[FALSE, ], phenotype_residual_SD = character()), "2 summary.parquet")
  error(aggregate(), "phenotype_residual_SD|numeric|schema")
  for (valid_scale in c(0, NA_real_)) {
    emit(sum = mutate(summary, phenotype_residual_SD = valid_scale))
    aggregate()
    stopifnot(identical(arrow::read_parquet("aggregate/susie.parquet")$phenotype_residual_SD, valid_scale))
  }
  emit()
  # IDs decoded from quoted BED fields must survive compact TSV aggregation.
  quote_id <- function(table) table %>%
    mutate(across(any_of(c("gene_id", "outcome_id", "phenotype_id", "molecular_trait_id")),
                  ~ replace(.x, .x == "a" | .x == "A", 'a"quoted')))
  readr::write_tsv(quote_id(wanted), "expected.tsv", quote = "needed")
  emit(idx = quote_id(index), sum = quote_id(summary),
       sets = quote_id(cs), diag = quote_id(diagnostic))
  aggregate()
  stopifnot(identical(read_tsv("aggregate/outcome_index.tsv")$phenotype_id,
                      c('a"quoted', "b", "c")),
            identical(read_tsv("aggregate/covariate_explained_variants.tsv")$molecular_trait_id,
                      'a"quoted'))
  readr::write_tsv(wanted, "expected.tsv")
  # Additional features can belong to this gene while being on another chromosome.
  for (variant_id in c("2_10_A_G", "chr2_10_A_G", "chr2_10_A_G_extra",
                       "chrchr2_10_A_G", "CHRchr2_10_A_G")) {
    emit(sum = mutate(summary, variant = variant_id, chromosome = "2"))
    aggregate()
    retained <- arrow::read_parquet("aggregate/susie.parquet")
    stopifnot(identical(retained$variant, variant_id), retained$chromosome == "2",
              retained$outcome_id == "a")
  }
  # Variant chromosome metadata must agree with the ID, not the gene chromosome.
  for (invalid_summary in list(
    mutate(summary, variant = "chr2_10_A_G"),
    mutate(summary, chromosome = NA_character_),
    mutate(summary, variant = NA_character_),
    mutate(summary, variant = "chr2"),
    mutate(summary, variant = "chr2_missing_A_G"),
    mutate(summary, variant = "chr2_10__G"),
    mutate(summary, variant = "chr2_10_A_"),
    mutate(summary, variant = "chr_10_A_G"))) {
    emit(sum = invalid_summary)
    error(aggregate(), "chromosome|variant|ID|autosome")
  }
  emit()
  # A complete status index cannot hide missing or unrelated high-purity variants.
  summary_cases <- list(
    missing_high_purity_summary = list(sum = summary[FALSE, ]),
    unknown_summary_set = list(sum = mutate(summary, cs_id = "unknown")),
    wrong_set_owner = list(sum = mutate(summary, outcome_id = "b", phenotype_id = "b",
                                        molecular_trait_id = "b")),
    summary_for_low_purity_set = list(sets = mutate(cs, low_purity = TRUE)),
    missing_set_member = list(sets = mutate(cs, cs_size = 2), sum = mutate(summary, cs_size = 2L)),
    wrong_summary_size = list(sum = mutate(summary, cs_size = 2L)),
    wrong_summary_flag = list(sum = mutate(summary, has_covariate_explained_variant = TRUE))
  )
  accepted_cases <- names(summary_cases)[vapply(summary_cases, function(arguments) {
    do.call(emit, arguments)
    result <- tryCatch({
      aggregate()
      NULL
    }, error = identity)
    !inherits(result, "error")
  }, logical(1))]
  if (length(accepted_cases)) stop("Aggregator accepted invalid credible-set summaries: ",
                                  paste(accepted_cases, collapse = ", "))
  # Low-purity sets remain in the index and CS table, with no summary variant rows.
  emit(sum = summary[FALSE, ], sets = mutate(cs, low_purity = TRUE))
  aggregate()
  stopifnot(nrow(arrow::read_parquet("aggregate/susie.parquet")) == 0L,
            nrow(arrow::read_parquet("aggregate/credible_sets.parquet")) == 1L)
  # A multi-variant set is valid only with all distinct variant rows present.
  two_variant_summary <- bind_rows(summary, mutate(summary, variant = "1_12_A_T", position = 12L)) %>%
    mutate(cs_size = 2L)
  emit(sum = two_variant_summary, sets = mutate(cs, cs_size = 2))
  aggregate()
  stopifnot(nrow(arrow::read_parquet("aggregate/susie.parquet")) == 2L,
    identical(arrow::read_parquet("aggregate/susie.parquet")$phenotype_residual_SD, c(1.25, 1.25)))
  for (bad_scale in c(2, NA_real_)) {
    inconsistent <- two_variant_summary
    inconsistent$phenotype_residual_SD[[2L]] <- bad_scale
    emit(sum = inconsistent, sets = mutate(cs, cs_size = 2))
    error(aggregate(), "phenotype_residual_SD|inconsistent")
  }
  emit()
  # Historical empty schemas have numeric ID/boolean columns; nonempty data stay numeric.
  empty_cs <- bind_cols(own[FALSE, ], InitEmptyCS())
  arrow::write_parquet(empty_cs, "2 cs.parquet")
  aggregate()
  stopifnot(is.numeric(arrow::read_parquet("aggregate/credible_sets.parquet")$cs_log10bf))
  emit(idx = index[-3, ])
  error(aggregate(), "missing|empty|batch")
  emit(idx = bind_rows(index, index[1, ]))
  error(aggregate(), "duplicate")
  emit(idx = mutate(index, status = "FAILED"))
  error(aggregate(), "FAILED")
  emit(sum = mutate(summary, phenotype_id = "foreign"))
  error(aggregate(), "ownership|foreign|expected")
  emit(sets = select(cs, -low_purity))
  error(aggregate(), "schema|column")
  emit(sets = mutate(cs, has_covariate_explained_variant = "FALSE"))
  error(aggregate(), "logical|boolean|type")
  emit(idx = mutate(index, converged = TRUE))
  error(aggregate(), "converged|status")
  emit(diag = mutate(diagnostic, gene_id = "foreign"))
  error(aggregate(), "ownership|foreign")
  emit(idx = mutate(index, chromosome = "chr2"))
  error(aggregate(), "chromosome")
  emit(sum = mutate(summary, variant = 123))
  error(aggregate(), "schema|text|ID")
  emit(idx = mutate(index, seed = 3e9))
  error(aggregate(), "seed")
  # All-empty summaries are valid only with complete, successful zero-CS statuses.
  zero_index <- mutate(index, credible_set_count = 0L)
  emit(idx = zero_index, sum = summary[FALSE, ], sets = empty_cs, diag = diagnostic[FALSE, ])
  aggregate()
  stopifnot(nrow(arrow::read_parquet("aggregate/susie.parquet")) == 0L,
            nrow(read_tsv("aggregate/outcome_index.tsv")) == 3L,
            is.double(arrow::read_parquet("aggregate/susie.parquet")$phenotype_residual_SD),
            is.logical(arrow::read_parquet("aggregate/susie.parquet")$feature_converged),
            is.logical(arrow::read_parquet("aggregate/credible_sets.parquet")$feature_converged))
  emit(idx = zero_index, sum = select(summary[FALSE, ], -phenotype_residual_SD),
       sets = empty_cs, diag = diagnostic[FALSE, ])
  aggregate()
  stopifnot(is.double(arrow::read_parquet("aggregate/susie.parquet")$phenotype_residual_SD))
  emit(idx = zero_index[-2, ], sum = summary[FALSE, ], sets = empty_cs, diag = diagnostic[FALSE, ])
  error(aggregate(), "Missing|missing")
  emit(idx = bind_rows(zero_index, zero_index[1, ]), sum = summary[FALSE, ],
       sets = empty_cs, diag = diagnostic[FALSE, ])
  error(aggregate(), "duplicate|Duplicate")
  emit()
  malformed <- readLines("1 index.tsv")
  malformed[[2]] <- paste0(malformed[[2]], "\textra-field")
  writeLines(malformed, "1 index.tsv")
  error(aggregate(), "parse|column|schema")
  emit()
  unlink("1 summary.parquet")
  error(aggregate(), "read|exist|local")
  error(aggregatePgenBatches("gs://bucket/expected.tsv", character(), character(),
                            character(), character(), "bad"), "localization")
  error(planPgenBatches("gs://bucket/phenotype.bed", "small.tsv", "chr1", 5,
                       "bad"), "localization")
  # The command entrypoints must pass named paths, including spaces, unchanged.
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
  writeLines(c("1 summary.parquet", "2 summary.parquet"), "summaries list.txt")
  writeLines(c("1 cs.parquet", "2 cs.parquet"), "sets list.txt")
  writeLines(c("1 index.tsv", "2 index.tsv"), "indices list.txt")
  writeLines(c("1 diagnostic.tsv", "2 diagnostic.tsv"), "diagnostics list.txt")
  run_cli("aggregate_pgen_susie.R", c("--expected-outcomes", "expected.tsv",
      "--summary-files", "summaries list.txt", "--credible-set-files", "sets list.txt",
      "--index-files", "indices list.txt", "--diagnostic-files", "diagnostics list.txt",
      "--output-dir", "cli aggregate"))
  stopifnot(nrow(read_tsv("cli aggregate/outcome_index.tsv")) == 3L)
  error(aggregatePgenBatches("expected.tsv", c("gs://bucket/a.parquet", "2 summary.parquet"),
                            c("1 cs.parquet", "2 cs.parquet"),
                            c("1 index.tsv", "2 index.tsv"),
                            c("1 diagnostic.tsv", "2 diagnostic.tsv"), "bad"), "localization")
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  run_pgen_batch_tests(normalizePath(if (length(args)) args[[1L]] else ".", mustWork = TRUE))
}
