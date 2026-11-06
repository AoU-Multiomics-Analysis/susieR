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
                     "SusieOutputFunctions.R", "PgenGeneFunctions.R")) {
    source(file.path(functions_path, filename), local = TRUE)
  }
  for (name in c("formatSusieOutputs", "writeSusieOutputs", "finemapPgenPhenotype",
                 "geneArtifactKey", "geneOutcomeSeed", "runPgenGene")) {
    if (!exists(name, mode = "function")) stop("Missing gene function: ", name)
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
    if (policy == "-1") gt[is.na(gt)] <- -1
    filtered <- filterMAF(gt, ancestry, MAF_threshold = maf, variant_list = variant_list)
    reference_input <- filtered
    if (tabix) {
      tsv <- file.path(temp, paste0("dosage-", policy, ".tsv"))
      dosage <- data.frame(CHROM = "1", POS = positions, REF = "A", ALT = "G") |>
        bind_cols(as.data.frame(if (policy == "-1") replace(raw, is.na(raw), -1) else raw))
      data.table::fwrite(dosage, tsv, sep = "\t", quote = FALSE, na = "NA", nThread = 1L)
      reference_input <- Rsamtools::bgzip(tsv, dest = paste0(tsv, ".bgz"), overwrite = TRUE)
      Rsamtools::indexTabix(reference_input, seq = 1L, start = 2L, end = 2L, skip = 1L, zeroBased = FALSE)
    }
    set.seed(181L)
    reference <- finemapPhenotype(trait, input_se, reference_input, cov, distance,
                                  ancestry_data = ancestry, MAF = maf, variant_list = variant_list,
                                  additional_genotypes = additional)
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
    result <- runPgenGene("gene/α", outcomes, reader, temp, 13L, "NA",
      provenance = list(inputs = list(pgen = list(basename = "fixture.pgen", sha256 = "supplied"))))
    stopifnot(identical(.Random.seed, rng), identical(RNGkind(), kind),
              identical(names(result), c("summary", "credible_sets", "diagnostics", "index", "files")),
              identical(names(result$index), c("gene_id", "gene_key", "outcome_id", "phenotype_id",
                "seed", "status", "converged", "credible_set_count", "reason")),
              is.integer(result$index$seed), is.integer(result$index$credible_set_count),
              is.logical(result$index$converged), all(file.exists(result$files)))
    fits <- readRDS(result$files[["rds"]])
    stopifnot(identical(names(fits), names(outcomes)),
              identical(unique(result$summary$molecular_trait_id), names(outcomes)[1:2]),
              result$index$credible_set_count[[3L]] == 0L,
              identical(unique(result$summary$region), unname(vapply(outcomes[1:2], `[[`, character(1), "region"))))
    position <- readr::read_tsv(result$files[["variant_positions"]], show_col_types = FALSE)
    stopifnot(identical(position$requested_window_start, c(-100, 120, -100)))
    provenance <- jsonlite::read_json(result$files[["provenance"]], simplifyVector = TRUE)
    stopifnot(identical(provenance$inputs$pgen$sha256, "supplied"),
              all(provenance$artifacts$basename %in% basename(result$files)))
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
    stopifnot(identical(tables$in_cs_variant_df, InitEmptyInCSVariantDf(TRUE)),
              identical(tables$lbf_df, InitEmptyLbfDf()),
              identical(tables$variant_df, InitEmptyVariantDf()),
              identical(tables$cs_df, InitEmptyCS()), nrow(tables$covariate_explained_variants) == 0L)
    files <- writeSusieOutputs(results, tables, file.path(temp, "empty"))
    stopifnot(identical(readRDS(files[["rds"]]), results), all(file.exists(files)))
    empty <- formatSusieOutputs(list(), tibble(phenotype_id = character(), region = character()))
    expect_identical_value(empty, tables, "empty result schemas")
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
