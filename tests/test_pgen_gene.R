# Real hardcalls and native fits test the per-gene interface.
run_pgen_gene_tests <- function(repo_dir) {
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(tidyr))
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(purrr))
  source(file.path(repo_dir, "tests", "test_helpers.R"), local = TRUE)
  functions_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH",
                               unset = file.path(repo_dir, "R", "utils"))
  for (filename in c("InitFunctions.R", "ImportFunctions.R", "MemoryHelpers.R",
                     "SusiePerformance.R", "SusieFunctions.R", "PgenFunctions.R",
                     "SusieOutputFunctions.R", "PgenGeneFunctions.R", "PgenBatchFunctions.R")) {
    source(file.path(functions_path, filename), local = TRUE)
  }
  for (name in c("formatSusieOutputs", "writeSusieOutputs", "finemapPgenPhenotype",
                 "geneArtifactKey", "geneOutcomeSeed", "runPgenGene")) {
    if (!exists(name, mode = "function")) stop("Missing gene function: ", name)
  }
  # Repeating extraction in the formatter must fail the call-count checks.
  # Changing output order, types, values, or serialization must fail parity.
  compare_reused_outputs <- function(results, regions, directory) {
    extracted <- purrr::map(results, extractResults)
    fresh <- formatSusieOutputs(results, regions)
    production_extract <- extractResults
    extraction_environment <- environment(formatSusieOutputs)
    extraction_calls <- 0L
    tryCatch({
      assign("extractResults", function(fitted) {
        extraction_calls <<- extraction_calls + 1L
        production_extract(fitted)
      }, envir = extraction_environment)
      reused <- formatSusieOutputs(results, regions, extracted_results = extracted)
    }, finally = {
      assign("extractResults", production_extract, envir = extraction_environment)
    })
    stopifnot(extraction_calls == 0L)
    expect_identical_value(reused, fresh, "reused output tables")
    expect_identical_value(
      formatSusieOutputs(results, regions, extracted_results = rev(extracted)),
      fresh, "named extraction order")
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    fresh_files <- writeSusieOutputs(results, fresh, file.path(directory, "fresh"))
    reused_files <- writeSusieOutputs(results, reused, file.path(directory, "reused"))
    expect_identical_value(names(reused_files), names(fresh_files), "output file types")
    expect_identical_value(readRDS(reused_files[["rds"]]), results, "complete saved fits")
    for (name in names(fresh_files)) {
      read_bytes <- function(path) readBin(path, "raw", n = file.info(path)$size)
      expect_identical_value(read_bytes(reused_files[[name]]),
                             read_bytes(fresh_files[[name]]), paste(name, "file bytes"))
    }
    for (name in c("summary", "lbf", "full")) {
      expect_identical_value(arrow::read_parquet(reused_files[[name]]),
                             arrow::read_parquet(fresh_files[[name]]),
                             paste(name, "decoded parquet"))
    }
    fresh
  }
  run_named_test("complete UTF-8 artifact keys and unambiguous stable seeds", {
    stopifnot(identical(geneArtifactKey("gene/é"), paste0("gene-", digest::digest(
      charToRaw(enc2utf8("gene/é")), algo = "sha256", serialize = FALSE))),
      nchar(geneArtifactKey("gene/é")) == 69L,
      identical(geneOutcomeSeed(1L, "a", "bc"), 161566456L),
      geneOutcomeSeed(1L, "a", "bc") != geneOutcomeSeed(1L, "ab", "c"),
      geneOutcomeSeed(1L, "a", "bc") != geneOutcomeSeed(2L, "a", "bc"))
    expect_error_message(geneArtifactKey(NA_character_), "gene_id")
    expect_error_message(geneOutcomeSeed(NA_integer_, "g", "o"), "base_seed")
  })
  temp <- tempfile("pgen-gene-tests-")
  dir.create(temp)
  on.exit(unlink(temp, recursive = TRUE), add = TRUE)
  run_named_test("failed gene status files preserve literal quotes in IDs", {
    quoted_gene <- 'g"quote'
    quoted_outcome <- 'outcome"quote'
    quoted_phenotype <- 'p"quote'
    quoted_se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = matrix(1, nrow = 1L, dimnames = list("other", "s1"))),
      rowData = data.frame(phenotype_id = "other", chromosome = "1", phenotype_pos = 1L))
    quoted_outcomes <- stats::setNames(list(list(
      phenotype_id = quoted_phenotype, se = quoted_se,
      covariates = matrix(numeric(), nrow = 1L, ncol = 0L),
      cis_distance = 5L, region = "1:1-6")), quoted_outcome)
    expect_error_message(runPgenGene(quoted_gene, quoted_outcomes, NULL, temp, 13L, "NA"),
                         "Expected one metadata row")
    prefix <- file.path(temp, "genes", geneArtifactKey(quoted_gene), geneArtifactKey(quoted_gene))
    index <- readr::read_tsv(paste0(prefix, ".outcome_index.tsv"),
      col_types = readr::cols(.default = readr::col_character()),
      na = character(), show_col_types = FALSE)
    expect_identical_value(index$gene_id, 'g"quote')
    expect_identical_value(index$outcome_id, 'outcome"quote')
    expect_identical_value(index$phenotype_id, 'p"quote')
    expect_identical_value(index$status, "FAILED")
    expect_identical_value(index$reason, 'Expected one metadata row for phenotype: p"quote')
  })
  run_named_test("gene diagnostic files preserve literal quotes in phenotype IDs", {
    tables <- formatSusieOutputs(list(), tibble::tibble(
      phenotype_id = character(), region = character()))
    tables$covariate_explained_variants <- tibble::tibble(
      molecular_trait_id = c('p"quote', '"pquote'), variant = c("1_10_A_G", "1_20_A_G"),
      input_scale = c(1, 2), residual_sd = c(1e-10, 4e-10),
      residual_sd_ratio = c(1e-10, 2e-10))
    files <- writeSusieOutputs(list(), tables, file.path(temp, "quoted-diagnostics"))
    diagnostics <- readr::read_tsv(files[["diagnostics"]],
      col_types = "ccddd", na = character(), show_col_types = FALSE)
    expect_identical_value(tibble::as_tibble(diagnostics), tables$covariate_explained_variants,
                           "diagnostic fields after TSV reading")
  })
  run_named_test("position summary appends preserve quotes and one header", {
    positions_file <- file.path(temp, "quoted-position-summary.tsv")
    first <- data.frame(phenotype_id = 'p"quote', context = "test", feature_position = 10,
      total_variants_tested = 2L, missing_text = NA_character_, stringsAsFactors = FALSE)
    second <- data.frame(phenotype_id = '"pquote', context = "test", feature_position = 20,
      total_variants_tested = 3L, missing_text = NA_character_, stringsAsFactors = FALSE)
    writeVariantPositionSummary(first, positions_file)
    writeVariantPositionSummary(second, positions_file)
    positions_read <- readr::read_tsv(positions_file, col_types = "ccdic",
      show_col_types = FALSE)
    expect_identical_value(positions_read$phenotype_id, c('p"quote', '"pquote'))
    expect_identical_value(positions_read$context, c("test", "test"))
    expect_identical_value(positions_read$feature_position, c(10, 20))
    expect_identical_value(positions_read$total_variants_tested, c(2L, 3L))
    expect_identical_value(positions_read$missing_text, c(NA_character_, NA_character_))
    expect_identical_value(sum(grepl('^"?phenotype_id"?\t', readLines(positions_file))), 1L)
    stopifnot(nrow(readr::problems(positions_read)) == 0L)
    ordinary <- first
    ordinary$phenotype_id <- "trait1"
    ordinary_file <- file.path(temp, "ordinary-position-summary.tsv")
    writeVariantPositionSummary(ordinary, ordinary_file)
    expect_identical_value(readLines(ordinary_file), c(
      "phenotype_id\tcontext\tfeature_position\ttotal_variants_tested\tmissing_text",
      "trait1\ttest\t10\t2\tNA"), "ordinary position summary bytes")
  })
  set.seed(913L)
  samples <- as.character(sample(seq_len(140L)))
  phenotype_samples <- samples[seq_len(100L)]
  positions <- seq_len(30L) * 10L
  ids <- paste0("chr1_", positions, "_A_G")
  raw <- matrix(as.double(rbinom(30L * 140L, 2L, 0.3)), nrow = 30L,
                dimnames = list(ids, samples))
  # Extra PSAM samples change ancestry MAF; first row also triggers dense math.
  raw[1L, ] <- rbinom(140L, 1L, 0.5)
  raw[2L, ] <- c(rep(0, 100L), rep(2, 40L))
  raw[2L, 1L] <- 1
  raw[3L, ] <- 0
  raw[3L, 2L] <- NA_real_
  raw[4L, ] <- NA_real_
  raw[5L, 3L] <- NA_real_
  paths <- file.path(temp, c("fixture.pgen", "fixture.pvar", "fixture.psam"))
  writer <- pgenlibr::NewPgenWriter(paths[[1L]], ncol(raw), nrow(raw), nonref_flags = FALSE)
  tryCatch({
    for (i in seq_len(nrow(raw))) pgenlibr::AppendBiallelic(writer, raw[i, ])
  }, finally = pgenlibr::ClosePgenWriter(writer))
  writeLines(c("#CHROM\tPOS\tID\tREF\tALT", paste("1", positions, ids, "A", "G", sep = "\t")), paths[[2L]])
  writeLines(c("#IID", samples), paths[[3L]])
  reader <- openPgenReader(paths[[1L]], paths[[2L]], paths[[3L]])
  on.exit(closePgenReader(reader), add = TRUE)
  covariates <- cbind(sex = raw[1L, phenotype_samples], pc = rnorm(100L))
  rownames(covariates) <- phenotype_samples
  y <- raw[8L, phenotype_samples] * 2 + rnorm(100L, sd = 0.3)
  y2 <- raw[15L, phenotype_samples] * 2 + rnorm(100L, sd = 0.3)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = rbind(trait1 = y, trait2 = y2)),
    colData = data.frame(sample_id = phenotype_samples, genotype_id = phenotype_samples,
                         row.names = phenotype_samples, qtl_group = "ALL"),
    rowData = data.frame(phenotype_id = c("trait1", "trait2"), chromosome = "1",
                         phenotype_pos = c(100L, 210L), row.names = c("trait1", "trait2")))
  covariates <- covariates[rev(phenotype_samples), , drop = FALSE]
  ancestry <- tibble(research_id = as.numeric(samples), ancestry_pred_other = "eur")
  # Compare the complete fitted object, including all native numerical fields.
  compare_path <- function(trait = "trait1", policy = "NA", maf = 0,
                           variant_list = NULL, additional = NULL, cov = covariates,
                           input_se = se, distance = 200L, tabix = FALSE) {
    meta <- as.data.frame(SummarizedExperiment::rowData(input_se))
    pos <- meta$phenotype_pos[meta$phenotype_id == trait]
    gt <- raw[positions >= max(1L, pos - distance) & positions <= pos + distance, , drop = FALSE]
    # Native imputation is independent of the legacy missing sentinel. Use the
    # NA reference so missing calls cannot cancel real ALT counts in filtering.
    reference_policy <- if (is.null(additional)) "NA" else policy
    if (reference_policy == "-1") gt[is.na(gt)] <- -1
    filtered <- filterMAF(gt, ancestry, MAF_threshold = maf, variant_list = variant_list)
    reference_input <- filtered
    if (tabix) {
      tsv <- file.path(temp, paste0("dosage-", policy, ".tsv"))
      dosage <- data.frame(CHROM = "1", POS = positions, REF = "A", ALT = "G") |>
        bind_cols(as.data.frame(if (reference_policy == "-1") replace(raw, is.na(raw), -1) else raw))
      data.table::fwrite(dosage, tsv, sep = "\t", quote = FALSE, na = "NA", nThread = 1L)
      reference_input <- Rsamtools::bgzip(tsv, dest = paste0(tsv, ".bgz"), overwrite = TRUE)
      Rsamtools::indexTabix(reference_input, seq = 1L, start = 2L, end = 2L, skip = 1L, zeroBased = FALSE)
    }
    # The PGEN path now discards all-missing variants under either missing-code
    # policy. Compare native and R imputation on that same retained variant set.
    reference_imputer <- if (is.null(additional)) {
      function(mat, missing_value) {
        keep <- rowSums(!is.na(mat) & mat != missing_value) > 0L
        impute_matrix_rowmean(mat[keep, , drop = FALSE], missing_value)
      }
    } else NULL
    set.seed(181L)
    reference <- finemapPhenotype(trait, input_se, reference_input, cov, distance,
                                  ancestry_data = ancestry, MAF = maf, variant_list = variant_list,
                                  additional_genotypes = additional,
                                  genotype_imputer = reference_imputer)
    expected_rng <- .Random.seed
    set.seed(181L)
    actual <- finemapPgenPhenotype(trait, input_se, reader, cov, distance, policy,
                                   ancestry_data = ancestry, MAF = maf, variant_list = variant_list,
                                   additional_genotypes = additional)
    expect_identical_value(actual, reference, "complete native fit")
    stopifnot(identical(expected_rng, .Random.seed))
    region <- tibble(phenotype_id = trait, region = paste0("chr1_", pos - distance, "_", pos + distance))
    expect_identical_value(formatSusieOutputs(setNames(list(actual), trait), region),
                           formatSusieOutputs(setNames(list(reference), trait), region), "decoded tables")
    actual
  }
  run_named_test("real PGEN matches raw matrix and Tabix under both missing policies", {
    for (policy in c("NA", "-1")) {
      invisible(compare_path(policy = policy))
      invisible(compare_path(policy = policy, tabix = TRUE, distance = 90L))
    }
  })
  run_named_test("PGEN fits impute on the first read without repeating R imputation", {
    reads <- new.env(parent = emptyenv())
    reads$native <- 0L
    reads$raw <- 0L
    reads$hardcall <- 0L
    reads$legacy <- 0L
    record <- function(field) substitute({
      .read_counter <- counter
      .read_counter[[key]] <- .read_counter[[key]] + 1L
    }, list(counter = reads, key = field))
    imputation_environment <- environment(finemapPhenotype)
    native_counter <- substitute({
      .read_counter <- counter
      .read_key <- if (meanimpute) "native" else "raw"
      .read_counter[[.read_key]] <- .read_counter[[.read_key]] + 1L
    }, list(counter = reads))
    suppressMessages(trace("ReadList", tracer = native_counter,
      where = asNamespace("pgenlibr"), print = FALSE))
    suppressMessages(trace("ReadIntList", tracer = record("hardcall"),
      where = asNamespace("pgenlibr"), print = FALSE))
    suppressMessages(trace("impute_matrix_rowmean", tracer = record("legacy"),
      where = imputation_environment, print = FALSE))
    tryCatch({
      # This window has no entirely missing variants, so no native retry is needed.
      invisible(compare_path(trait = "trait2", distance = 70L))
      stopifnot(reads$native == 1L, reads$hardcall == 0L,
                reads$raw == 1L,
                reads$legacy == 1L) # R imputation runs only in the reference fit.
    }, finally = {
      suppressMessages(untrace("ReadList", where = asNamespace("pgenlibr")))
      suppressMessages(untrace("ReadIntList", where = asNamespace("pgenlibr")))
      suppressMessages(untrace("impute_matrix_rowmean", where = imputation_environment))
    })
  })
  run_named_test("full-cohort ancestry MAF and external list precedence and order", {
    full_maf <- filterMAF(raw, ancestry, MAF_threshold = 0.1)
    subset_maf <- filterMAF(raw[, phenotype_samples], ancestry, MAF_threshold = 0.1)
    stopifnot(ids[[2L]] %in% rownames(full_maf), !ids[[2L]] %in% rownames(subset_maf))
    invisible(compare_path(maf = 0.1))
    allowlist <- file.path(temp, "variants.tsv")
    readr::write_tsv(tibble(variant = rev(ids[c(1L, 8L, 15L, 20L)])), allowlist)
    fit <- compare_path(maf = 0.49, variant_list = allowlist)
    stopifnot(identical(fit$variant_id, ids[c(1L, 8L, 15L, 20L)]))
    # The legacy filter drops dimensions when exactly one variant survives.
    # Preserve its error until the legacy path is changed in a separate task.
    readr::write_tsv(tibble(variant = ids[[8L]]), allowlist)
    one_variant <- filterMAF(raw, ancestry, variant_list = allowlist)
    reference_error <- tryCatch(
      finemapPhenotype("trait1", se, one_variant, covariates, 200L),
      error = conditionMessage)
    adapter_error <- tryCatch(
      finemapPgenPhenotype("trait1", se, reader, covariates, 200L, "NA",
                           variant_list = allowlist), error = conditionMessage)
    stopifnot(is.character(reference_error), identical(adapter_error, reference_error))
    subset_reader <- openPgenReader(paths[[1L]], paths[[2L]], paths[[3L]], phenotype_samples)
    tryCatch(expect_error_message(finemapPgenPhenotype("trait1", se, subset_reader,
      covariates, 200L, "NA", ancestry, 0.1), "complete"),
      finally = closePgenReader(subset_reader))
  })
  run_named_test("native-first counts retain physical duplicate labels through filtering", {
    duplicate_positions <- positions
    duplicate_positions[[5L]] <- positions[[1L]]
    duplicate_pvar <- file.path(temp, "duplicate-counts.pvar")
    writeLines(c("#CHROM\tPOS\tID\tREF\tALT",
      paste("1", duplicate_positions, ids, "A", "G", sep = "\t")), duplicate_pvar)
    duplicate_reader <- openPgenReader(paths[[1L]], duplicate_pvar, paths[[3L]])
    tryCatch({
      fitted <- finemapPgenPhenotype("trait1", se, duplicate_reader, covariates, 200L, "NA")
      duplicate_rows <- fitted$variant_id == ids[[1L]]
      expect_identical_value(fitted$cohort_an[duplicate_rows], c(200, 198),
                             "duplicate rows keep different observed sample counts")
      expected_alt <- rowSums(raw[c(1L, 5L), phenotype_samples, drop = FALSE], na.rm = TRUE)
      expected_mac <- unname(pmin(expected_alt, c(200, 198) - expected_alt))
      expect_identical_value(fitted$cohort_mac[duplicate_rows], expected_mac,
                             "duplicate rows keep different observed allele counts")
      allowlist <- file.path(temp, "duplicate-counts-allowlist.tsv")
      readr::write_tsv(tibble(variant = ids[c(1L, 8L)]), allowlist)
      listed <- finemapPgenPhenotype("trait1", se, duplicate_reader, covariates, 200L,
                                     "NA", variant_list = allowlist)
      expect_identical_value(listed$variant_id, ids[c(1L, 8L)], "existing first-match list order")
      expect_identical_value(listed$cohort_an, c(200, 200), "list uses the matching physical rows")
    }, finally = closePgenReader(duplicate_reader))
  })
  run_named_test("additional samples and genotype/response dense fallback retain every field", {
    additional_samples <- phenotype_samples[seq_len(80L)]
    additional <- matrix(as.double(rbinom(2L * 80L, 2L, 0.3)), nrow = 2L,
      dimnames = list(c("chr1_310_A_T", "chr1_320_A_T"), rev(additional_samples)))
    invisible(compare_path(additional = additional))
    fit <- compare_path()
    stopifnot(ids[[1L]] %in% fit$covariate_explained_variants$variant_id)
    transformed_y <- qnorm((rank(y) - 0.5) / length(y))
    near_cov <- cbind(pc = transformed_y + 1e-10 * sin(seq_along(y)))
    rownames(near_cov) <- phenotype_samples
    invisible(compare_path(cov = near_cov))
    na_se <- se
    SummarizedExperiment::assay(na_se)[1L, 2L] <- NA_real_
    set.seed(181L)
    raw_error <- tryCatch(finemapPhenotype("trait1", na_se, raw, covariates, 200L),
                           error = conditionMessage)
    set.seed(181L)
    pgen_error <- tryCatch(finemapPgenPhenotype("trait1", na_se, reader, covariates, 200L, "NA"),
                            error = conditionMessage)
    expect_identical_value(pgen_error, raw_error, "NA phenotype behavior")
  })
  run_named_test("timing callbacks cover loading preparation and fitting without changing fits or RNG", {
    # Removing a stage callback or storing measurements in a fit must fail this check.
    events <- list()
    record_stage <- function(stage, event) {
      events[[length(events) + 1L]] <<- c(stage = stage, event = event)
    }
    set.seed(181L)
    reference <- finemapPgenPhenotype("trait1", se, reader, covariates, 200L, "NA")
    expected_rng <- .Random.seed
    set.seed(181L)
    timed <- finemapPgenPhenotype("trait1", se, reader, covariates, 200L, "NA",
                                  timing_callback = record_stage)
    expect_identical_value(timed, reference, "timing leaves complete fit unchanged")
    expect_identical_value(.Random.seed, expected_rng, "timing leaves RNG unchanged")
    expect_identical_value(events, list(
      c(stage = "data_loading", event = "start"),
      c(stage = "data_loading", event = "complete"),
      c(stage = "data_preparation", event = "start"),
      c(stage = "data_preparation", event = "complete"),
      c(stage = "fine_mapping", event = "start"),
      c(stage = "fine_mapping", event = "complete")))
    events <- list()
    invalid_covariates <- cbind(covariates, duplicate_pc = covariates[, "pc"])
    expect_error_message(finemapPgenPhenotype("trait1", se, reader,
      invalid_covariates, 200L, "NA", timing_callback = record_stage), "singular")
    stopifnot(identical(tail(events, 1L)[[1L]],
                       c(stage = "data_preparation", event = "failed")),
              !any(vapply(events, function(event) event[["stage"]] == "fine_mapping", logical(1))))
    reusable <- prepareReusableGenotype("trait1", se, raw, covariates, 200L)
    events <- list()
    set.seed(181L)
    reference <- finemapPhenotype("trait1", se, raw, covariates, 200L,
                                  reusable_genotype = reusable)
    expected_rng <- .Random.seed
    set.seed(181L)
    timed <- finemapPhenotype("trait1", se, raw, covariates, 200L,
                               reusable_genotype = reusable, timing_callback = record_stage)
    expect_identical_value(timed, reference, "timing leaves reusable fit unchanged")
    expect_identical_value(.Random.seed, expected_rng, "reusable timing leaves RNG unchanged")
    expect_identical_value(events, list(
      c(stage = "data_preparation", event = "start"),
      c(stage = "data_preparation", event = "complete"),
      c(stage = "fine_mapping", event = "start"),
      c(stage = "fine_mapping", event = "complete")))
  })
  outcomes <- list(
    `outcome / α` = list(phenotype_id = "trait1", se = se, covariates = covariates,
                         cis_distance = 200L, region = "chr1_-100_300"),
    outcome2 = list(phenotype_id = "trait2", se = se, covariates = covariates,
                    cis_distance = 90L, region = "chr1_120_300"),
    no_sets = list(phenotype_id = "trait1", se = se, covariates = near_cov,
                   cis_distance = 200L, region = "chr1_-100_300"))
  run_named_test("multi-outcome artifacts, deterministic order, compact return and RNG restore", {
    set.seed(997L)
    rng <- .Random.seed
    kind <- RNGkind()
    production_extract <- extractResults
    extracted_fits <- list()
    tryCatch({
      extractResults <- function(fitted) {
        extracted_fits[[length(extracted_fits) + 1L]] <<- fitted
        production_extract(fitted)
      }
      result <- runPgenGene("gene/α", outcomes, reader, temp, 13L, "NA",
        provenance = list(inputs = list(pgen = list(basename = "fixture.pgen", sha256 = "supplied"))))
    }, finally = {
      extractResults <- production_extract
    })
    if (length(extracted_fits) != length(outcomes)) {
      stop("Expected one extractResults() call per outcome; expected ", length(outcomes),
           " calls, received ", length(extracted_fits), call. = FALSE)
    }
    stopifnot(identical(.Random.seed, rng), identical(RNGkind(), kind),
              identical(names(result), c("summary", "credible_sets", "diagnostics", "index", "files")),
              identical(names(result$index), c("gene_id", "gene_key", "outcome_id", "phenotype_id",
                "seed", "status", "converged", "credible_set_count", "reason",
                "data_loading_seconds", "data_preparation_seconds", "fine_mapping_seconds",
                "outcome_total_seconds")),
              is.integer(result$index$seed), is.integer(result$index$credible_set_count),
              is.logical(result$index$converged), all(file.exists(result$files)))
    timing_columns <- c("data_loading_seconds", "data_preparation_seconds", "fine_mapping_seconds",
                        "outcome_total_seconds")
    stopifnot(all(vapply(result$index[timing_columns], is.double, logical(1))),
              all(is.finite(as.matrix(result$index[timing_columns]))),
              all(as.matrix(result$index[timing_columns]) >= 0),
              all(rowSums(as.matrix(result$index[timing_columns[1:3]])) <=
                    result$index$outcome_total_seconds))
    fits <- readRDS(result$files[["rds"]])
    expect_identical_value(extracted_fits, unname(fits), "one extraction for each saved fit")
    stopifnot(identical(names(fits), names(outcomes)),
              identical(unique(result$summary$molecular_trait_id), names(outcomes)[1:2]),
              result$index$credible_set_count[[3L]] == 0L,
              identical(unique(result$summary$region), unname(vapply(outcomes[1:2], `[[`, character(1), "region"))))
    regions <- tibble(phenotype_id = names(outcomes),
                      region = vapply(outcomes, `[[`, character(1), "region"))
    fresh <- compare_reused_outputs(fits, regions, file.path(temp, "reuse-mixed"))
    expect_identical_value(result$summary, fresh$in_cs_variant_df, "gene summary rows")
    expect_identical_value(result$credible_sets, fresh$cs_df, "gene credible-set rows")
    expect_identical_value(result$diagnostics, fresh$covariate_explained_variants,
                           "gene diagnostic rows")
    counts <- vapply(fits, function(fitted) {
      sets <- extractResults(fitted)$cs_df
      if (is.null(sets)) 0L else as.integer(nrow(sets))
    }, integer(1))
    expect_identical_value(result$index$credible_set_count, unname(counts), "status set counts")
    extracted <- purrr::map(fits, extractResults)
    expect_error_message(formatSusieOutputs(fits, regions,
      extracted_results = unname(extracted)), "named list")
    expect_error_message(formatSusieOutputs(fits, regions,
      extracted_results = extracted[-1L]), "each outcome")
    extra <- extracted
    extra$unrequested <- extracted[[1L]]
    expect_error_message(formatSusieOutputs(fits, regions,
      extracted_results = extra), "each outcome")
    duplicate_extracted <- extracted
    names(duplicate_extracted)[[2L]] <- names(duplicate_extracted)[[1L]]
    expect_error_message(formatSusieOutputs(fits, regions,
      extracted_results = duplicate_extracted), "unique")
    incomplete_extracted <- extracted
    incomplete_extracted[[1L]]$lbf_df <- NULL
    expect_error_message(formatSusieOutputs(fits, regions,
      extracted_results = incomplete_extracted), "cs_df, variant_df, and lbf_df")
    invalid_extracted <- extracted
    invalid_extracted[[1L]]$cs_df <- "invalid table"
    expect_error_message(formatSusieOutputs(fits, regions,
      extracted_results = invalid_extracted), "data frames or NULL")
    all_flagged <- purrr::map(fits, function(fitted) {
      fitted$covariate_explained_variants <- tibble(
        variant_id = fitted$variant_id, input_scale = 1,
        residual_sd = 1e-10, residual_sd_ratio = 1e-10)
      fitted
    })
    flagged_tables <- compare_reused_outputs(all_flagged, regions, file.path(temp, "reuse-all-flagged"))
    stopifnot(nrow(flagged_tables$cs_df) > 0L,
              all(flagged_tables$cs_df$has_covariate_explained_variant),
              all(flagged_tables$in_cs_variant_df$has_covariate_explained_variant))
    old_fits <- purrr::map(fits, function(fitted) {
      fitted$covariate_explained_variants <- NULL
      fitted
    })
    old_tables <- compare_reused_outputs(old_fits, regions, file.path(temp, "reuse-old-fits"))
    stopifnot(!any(old_tables$cs_df$has_covariate_explained_variant),
              nrow(old_tables$covariate_explained_variants) == 0L)
    position <- readr::read_tsv(result$files[["variant_positions"]], show_col_types = FALSE)
    stopifnot(identical(position$requested_window_start, c(-100, 120, -100)))
    provenance <- jsonlite::read_json(result$files[["provenance"]], simplifyVector = TRUE)
    stopifnot(identical(provenance$inputs$pgen$sha256, "supplied"),
              all(provenance$artifacts$basename %in% basename(result$files)))
    for (column in timing_columns[1:3]) {
      stopifnot(abs(provenance$timing[[column]] - sum(result$index[[column]])) < 1e-6)
    }
    stopifnot(is.finite(provenance$timing$gene_total_seconds),
              provenance$timing$gene_total_seconds >= sum(result$index$outcome_total_seconds),
              identical(provenance$timing$scope,
                        "gene processing and output writing; excludes final hashing and provenance writing"),
              any(grepl("Gene timing:", readLines(result$files[["log"]]), fixed = TRUE)))
    # A log update after artifact hashing would make the published hash stale.
    log_record <- provenance$artifacts[provenance$artifacts$name == "log", ]
    stopifnot(identical(log_record$sha256,
                       digest::digest(file = result$files[["log"]], algo = "sha256")))
    # The aggregator must parse timing fields, preserve their values, and reject text in a numeric field.
    aggregate_root <- file.path(temp, "timing-aggregate")
    dir.create(aggregate_root)
    batch_index <- result$index[1:2, ] |>
      dplyr::mutate(batch_id = "batch_000001", chromosome = "1")
    expected <- batch_index |>
      dplyr::select(batch_id, chromosome, gene_id, outcome_id, phenotype_id) |>
      dplyr::mutate(phenotype_pos = c(100, 210), requested_window_start = c(-100, 120),
                    requested_window_end = c(300, 300), window_start = c(1, 120),
                    window_end = 300, region = c("chr1_-100_300", "chr1_120_300"))
    mapping <- dplyr::select(expected, batch_id, gene_id, outcome_id, phenotype_id)
    decorate <- function(table) table |>
      dplyr::filter(.data$molecular_trait_id %in% mapping$outcome_id) |>
      dplyr::left_join(mapping, by = c("molecular_trait_id" = "outcome_id")) |>
      dplyr::mutate(outcome_id = molecular_trait_id)
    aggregate_inputs <- file.path(aggregate_root, c("expected.tsv", "summary.parquet",
      "cs.parquet", "index.tsv", "diagnostics.tsv"))
    readr::write_tsv(expected, aggregate_inputs[[1L]])
    arrow::write_parquet(decorate(result$summary), aggregate_inputs[[2L]])
    arrow::write_parquet(decorate(result$credible_sets), aggregate_inputs[[3L]])
    readr::write_tsv(batch_index, aggregate_inputs[[4L]])
    diagnostics <- decorate(result$diagnostics) |>
      dplyr::mutate(chromosome = "1")
    readr::write_tsv(diagnostics, aggregate_inputs[[5L]])
    aggregate_run <- function() aggregatePgenBatches(aggregate_inputs[[1L]],
      aggregate_inputs[[2L]], aggregate_inputs[[3L]], aggregate_inputs[[4L]],
      aggregate_inputs[[5L]], file.path(aggregate_root, "combined"))
    aggregate_files <- aggregate_run()
    combined <- readr::read_tsv(aggregate_files[["index"]], show_col_types = FALSE)
    for (column in timing_columns) {
      stopifnot(is.double(combined[[column]]),
                all(abs(combined[[column]] - batch_index[[column]]) < 1e-9))
    }
    invalid_index <- batch_index
    invalid_index$fine_mapping_seconds <- c("invalid", "1")
    readr::write_tsv(invalid_index, aggregate_inputs[[4L]])
    expect_error_message(aggregate_run(), "Invalid numeric column: fine_mapping_seconds")
    negative_index <- batch_index
    negative_index$data_loading_seconds[[1L]] <- -1
    readr::write_tsv(negative_index, aggregate_inputs[[4L]])
    negative_files <- aggregate_run()
    negative_combined <- readr::read_tsv(negative_files[["index"]], show_col_types = FALSE)
    stopifnot(negative_combined$data_loading_seconds[[1L]] == -1,
              identical(negative_combined$data_loading_seconds,
                        negative_index$data_loading_seconds))
    reversed <- runPgenGene("gene/α", rev(outcomes), reader, file.path(temp, "reverse"), 13L, "NA")
    reverse_fits <- readRDS(reversed$files[["rds"]])
    for (id in names(outcomes)) expect_identical_value(fits[[id]], reverse_fits[[id]], "order-independent fit")
    single <- runPgenGene("gene/α", outcomes[1L], reader, file.path(temp, "single"), 13L, "NA")
    expect_identical_value(fits[[1L]], readRDS(single$files[["rds"]])[[1L]], "single fit")
    # A caller without RNG state must remain without RNG state.
    rm(".Random.seed", envir = .GlobalEnv)
    invisible(runPgenGene("other gene", outcomes[1L], reader, temp, 13L, "NA"))
    stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    assign(".Random.seed", rng, envir = .GlobalEnv)
    duplicate <- outcomes
    names(duplicate) <- rep("duplicate", 2L)
    expect_error_message(runPgenGene("bad", duplicate, reader, temp, 1L, "NA"), "unique")
  })
  run_named_test("zero-CS schemas and old fits remain exportable without quitting", {
    set.seed(191L)
    fit <- finemapPgenPhenotype("trait1", se, reader, covariates, 200L, "NA")
    fit$sets$cs <- NULL
    fit$sets$cs_index <- integer()
    fit$sets$purity <- fit$sets$purity[0, , drop = FALSE]
    fit$covariate_explained_variants <- NULL
    results <- list(old = fit)
    tables <- formatSusieOutputs(results, tibble(phenotype_id = "old", region = "chr1_-100_300"))
    stopifnot(identical(tables$in_cs_variant_df,
                        mutate(InitEmptyInCSVariantDf(TRUE, include_phenotype_residual_sd = TRUE,
                                                     include_cohort_maf = TRUE,
                                                     include_cohort_counts = TRUE),
                               feature_converged = logical())),
              identical(tables$lbf_df, InitEmptyLbfDf()),
              identical(tables$variant_df,
                        mutate(InitEmptyVariantDf(include_phenotype_residual_sd = TRUE),
                               feature_converged = logical())),
              identical(tables$cs_df, mutate(InitEmptyCS(), feature_converged = logical())),
              nrow(tables$covariate_explained_variants) == 0L)
    files <- writeSusieOutputs(results, tables, file.path(temp, "empty"))
    stopifnot(identical(readRDS(files[["rds"]]), results), all(file.exists(files)))
    empty <- formatSusieOutputs(list(), tibble(phenotype_id = character(), region = character()))
    expect_identical_value(empty, tables, "empty result schemas")
    expect_identical_value(compare_reused_outputs(results,
      tibble(phenotype_id = "old", region = "chr1_-100_300"),
      file.path(temp, "reuse-no-sets")), tables, "reused no-set schemas")
    expect_identical_value(compare_reused_outputs(list(),
      tibble(phenotype_id = character(), region = character()),
      file.path(temp, "reuse-empty")), empty, "reused empty schemas")
    expect_error_message(formatSusieOutputs(setNames(list(fit, fit), c("x", "x")),
      tibble(phenotype_id = "x", region = "r")), "unique")
  })
  run_named_test("failed fits write status and headers without a success manifest", {
    invalid <- outcomes[1L]
    invalid[[1L]]$phenotype_id <- "absent"
    root <- file.path(temp, "failure")
    expect_error_message(runPgenGene("failed", invalid, reader, root, 13L, "NA"), "absent")
    prefix <- file.path(root, "genes", geneArtifactKey("failed"), geneArtifactKey("failed"))
    index <- readr::read_tsv(paste0(prefix, ".outcome_index.tsv"), show_col_types = FALSE)
    stopifnot(identical(index$status, "FAILED"), !file.exists(paste0(prefix, ".provenance.json")),
              length(readLines(paste0(prefix, ".variant_position_summary.tsv"))) == 1L,
              length(readLines(paste0(prefix, ".covariate_explained_variants.tsv"))) == 1L,
              any(grepl("FAILED", readLines(paste0(prefix, ".log")))))
    stopifnot(is.finite(index$data_loading_seconds), index$data_loading_seconds >= 0,
              is.na(index$data_preparation_seconds), is.na(index$fine_mapping_seconds),
              is.finite(index$outcome_total_seconds), index$outcome_total_seconds >= index$data_loading_seconds,
              any(grepl("Gene timing:", readLines(paste0(prefix, ".log")), fixed = TRUE)))
    constant <- outcomes[1L]
    constant_se <- se
    SummarizedExperiment::assay(constant_se)[1L, ] <- 1
    constant[[1L]]$se <- constant_se
    expect_error_message(runPgenGene("constant", constant, reader, root, 13L, "NA"), "failed:")
    prefix <- file.path(root, "genes", geneArtifactKey("constant"), geneArtifactKey("constant"))
    index <- readr::read_tsv(paste0(prefix, ".outcome_index.tsv"), show_col_types = FALSE)
    stages <- c("data_loading_seconds", "data_preparation_seconds", "fine_mapping_seconds")
    stopifnot(identical(index$status, "FAILED"),
              all(is.finite(as.matrix(index[stages]))), all(as.matrix(index[stages]) >= 0),
              sum(as.matrix(index[stages])) <= index$outcome_total_seconds,
              any(grepl("stage: fine_mapping; status: failed", readLines(paste0(prefix, ".log")), fixed = TRUE)))
  })
  run_named_test("NONCONVERGED native fits stay in the index and RDS", {
    # A one-iteration native fit makes status persistence deterministic. The
    # full adapter and fixed fitting settings are checked above without a stub.
    set.seed(37L)
    status_x <- matrix(rnorm(60L * 12L), nrow = 60L,
      dimnames = list(NULL, paste0("chr1_", seq_len(12L), "_A_G")))
    incomplete <- suppressWarnings(susieR::susie(status_x,
      cbind(2 * status_x[, 2L] + rnorm(60L)), L = 10, max_iter = 1,
      compute_univariate_zscore = TRUE))
    incomplete$variant_id <- colnames(status_x)
    stopifnot(identical(incomplete$converged, FALSE))
    production_adapter <- finemapPgenPhenotype
    tryCatch({
      finemapPgenPhenotype <- function(...) incomplete
      result <- runPgenGene("incomplete", outcomes[1L], reader, temp, 13L, "NA")
      stopifnot(identical(result$index$status, "NONCONVERGED"),
                identical(result$index$converged, FALSE),
                all(is.na(result$index[c("data_loading_seconds", "data_preparation_seconds", "fine_mapping_seconds")])),
                is.finite(result$index$outcome_total_seconds),
                identical(readRDS(result$files[["rds"]])[[1L]], incomplete))
    }, finally = {
      finemapPgenPhenotype <- production_adapter
    })
  })
  run_named_test("RNG kind and state restore after success and error", {
    caller_kind <- RNGkind()
    caller_state <- .Random.seed
    tryCatch({
      RNGkind("L'Ecuyer-CMRG")
      set.seed(51L)
      expected_kind <- RNGkind()
      expected_state <- .Random.seed
      invisible(.pgenGeneWithSeed(13L, rnorm(4L)))
      stopifnot(identical(RNGkind(), expected_kind), identical(.Random.seed, expected_state))
      expect_error_message(.pgenGeneWithSeed(13L, stop("fit error")), "fit error")
      stopifnot(identical(RNGkind(), expected_kind), identical(.Random.seed, expected_state))
    }, finally = {
      do.call(RNGkind, as.list(caller_kind))
      assign(".Random.seed", caller_state, envir = .GlobalEnv)
    })
  })
  source(file.path(repo_dir, "tests", "test_covariate_explained_outputs.R"), local = TRUE)
  run_covariate_explained_output_tests(repo_dir)
  message("All PGEN gene tests passed.")
  invisible(TRUE)
}
