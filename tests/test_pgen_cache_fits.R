# Real fits check that prepared caches preserve the scientific results, output
# tables, and RNG state. Window order includes overlap, a gap, and a return.
run_pgen_cache_fit_tests <- function(repo_dir) {
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
  cache_source <- file.path(functions_path, "PgenCacheFunctions.R")
  if (file.exists(cache_source)) source(cache_source, local = TRUE)

  directory <- tempfile("pgen-cache-fits-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  set.seed(184L)
  samples <- as.character(sample(1001:1100))
  selected_samples <- samples[rev(seq_len(80L))]
  positions <- seq_len(40L) * 100L
  variant_ids <- paste0("chr1_", positions, "_A_G")
  genotypes <- matrix(as.double(rbinom(40L * length(samples), 2L, 0.3)),
                      nrow = 40L, dimnames = list(variant_ids, samples))
  genotypes[3L, match(selected_samples[[4L]], samples)] <- NA_real_
  genotypes[9L, match(selected_samples[[8L]], samples)] <- NA_real_
  genotypes[31L, ] <- NA_real_
  genotypes[39L, ] <- 0
  files <- file.path(directory, c("input.pgen", "input.pvar", "input.psam"))
  writer <- pgenlibr::NewPgenWriter(files[[1L]], length(samples), nrow(genotypes),
                                   nonref_flags = FALSE)
  tryCatch({
    for (i in seq_len(nrow(genotypes))) {
      pgenlibr::AppendBiallelic(writer, genotypes[i, ])
    }
  }, finally = pgenlibr::ClosePgenWriter(writer))
  writeLines(c("#CHROM\tPOS\tID\tREF\tALT",
    paste("1", positions, variant_ids, "A", "G", sep = "\t")), files[[2L]])
  writeLines(c("#IID", samples), files[[3L]])
  reader <- do.call(openPgenReader, as.list(files))
  on.exit(closePgenReader(reader), add = TRUE)

  phenotype_ids <- c("left", "overlap", "distant", "return")
  feature_positions <- c(700L, 1000L, 3500L, 750L)
  causal_rows <- c(7L, 10L, 35L, 7L)
  phenotype_values <- t(vapply(causal_rows, function(i) {
    2 * genotypes[i, selected_samples] + rnorm(length(selected_samples), sd = 0.4)
  }, numeric(length(selected_samples))))
  dimnames(phenotype_values) <- list(phenotype_ids, selected_samples)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = phenotype_values),
    colData = data.frame(sample_id = selected_samples, genotype_id = selected_samples,
      qtl_group = "ALL", row.names = selected_samples),
    rowData = data.frame(phenotype_id = phenotype_ids, chromosome = "1",
      phenotype_pos = feature_positions, row.names = phenotype_ids))
  covariates <- cbind(pc1 = rnorm(length(selected_samples)),
                       pc2 = rnorm(length(selected_samples)))
  rownames(covariates) <- selected_samples
  # Supplying covariates in another order must not change the aligned design.
  covariates <- covariates[rev(selected_samples), , drop = FALSE]
  distance <- 500L
  worklist <- tibble::tibble(batch_id = "batch_000001", chromosome = "1",
    gene_id = phenotype_ids, outcome_id = phenotype_ids, phenotype_id = phenotype_ids,
    phenotype_pos = feature_positions,
    requested_window_start = feature_positions - distance,
    requested_window_end = feature_positions + distance,
    window_start = pmax(1L, feature_positions - distance),
    window_end = feature_positions + distance,
    region = paste0("chr1:", pmax(1L, feature_positions - distance), "-",
                    feature_positions + distance))

  observed_inputs <- NULL
  original_fit <- fitSusieModel
  fitSusieModel <- function(X, y, ...) {
    observed_inputs <<- list(X = X, y = y)
    original_fit(X, y, ...)
  }
  fit <- function(id, cov = covariates, cache = NULL, explicit = TRUE,
                   timing_callback = NULL) {
    set.seed(259L)
    arguments <- list(phenotype_id = id, se = se, reader = reader, covariates = cov,
                      cis_distance = distance, missing_code = "NA")
    if (explicit) {
      arguments["genotype_cache"] <- list(cache)
      arguments$cache_outcome_id <- id
    }
    if (!is.null(timing_callback)) arguments$timing_callback <- timing_callback
    result <- do.call(finemapPgenPhenotype, arguments)
    list(result = result, rng = .Random.seed, inputs = observed_inputs,
      outcome_id = id, cache_mode = if (is.null(cache)) "none" else cache$mode)
  }
  make_cache <- function(mode, cov = covariates) {
    plan <- planPgenGenotypeCache(reader, worklist, selected_samples, mode = mode,
      cache_memory_gb = 0.1, job_memory_gb = 1)
    createPgenGenotypeCache(reader, plan, selected_samples,
      cbind(cov[selected_samples, , drop = FALSE], 1))
  }
  # Cache block widths can change BLAS rounding. Native Brent prior-variance
  # optimization stops at sqrt(machine epsilon), so complete fitted parameters
  # can differ more than the inputs. These are limits for this fixture, not a
  # bound for arbitrary data. Guarded fits still require exact equality.
  model_tolerance <- 1e-7
  pip_absolute_tolerance <- 1e-8
  compare_model <- function(actual, reference, exact = FALSE) {
    expect_identical_value(actual$variant_id, reference$variant_id, "fit variant order")
    expect_identical_value(actual$sets$cs, reference$sets$cs, "credible-set memberships")
    expect_identical_value(actual$converged, reference$converged, "convergence status")
    expect_identical_value(actual$niter, reference$niter, "iteration count")
    expect_identical_value(length(actual$pip), length(reference$pip), "PIP count")
    difference <- abs(actual$pip - reference$pip)
    if (any(!is.finite(difference)) || max(difference) > pip_absolute_tolerance) {
      stop("PIP absolute difference exceeds 1e-8: ",
        format(max(difference), digits = 16L), call. = FALSE)
    }
    if (exact) {
      expect_identical_value(actual, reference, "guarded complete fit")
    } else {
      comparison <- all.equal(actual, reference, tolerance = model_tolerance,
                              check.attributes = TRUE)
      if (!isTRUE(comparison)) stop(paste(comparison, collapse = "\n"), call. = FALSE)
    }
  }
  compare_fit <- function(actual, reference, exact = FALSE) {
    input_differences <- vapply(c("X", "y"), function(name) {
      left <- actual$inputs[[name]]
      right <- reference$inputs[[name]]
      expect_identical_value(attributes(left), attributes(right),
        paste("Adjusted fit input", name, "attributes"))
      difference <- abs(left - right)
      if (any(!is.finite(difference)) || max(difference) > 1e-13) {
        stop("Adjusted fit input ", name, " differs for ", actual$outcome_id,
          "; max absolute difference=", format(max(difference), digits = 16L), call. = FALSE)
      }
      max(difference)
    }, numeric(1))
    if (!exact && actual$cache_mode != "none") {
      message("Cache fit diagnostics: outcome=", actual$outcome_id,
        "; mode=", actual$cache_mode,
        "; X_max_abs=", format(input_differences[["X"]], digits = 16L),
        "; y_max_abs=", format(input_differences[["y"]], digits = 16L),
        "; PIP_max_abs=", format(max(abs(actual$result$pip - reference$result$pip)), digits = 16L),
        "; mu_max_abs=", format(max(abs(actual$result$mu - reference$result$mu)), digits = 16L))
    }
    expect_identical_value(actual$rng, reference$rng, "fit RNG state")
    compare_model(actual$result, reference$result, exact)
  }
  contains_environment <- function(value) {
    is.environment(value) ||
      (is.list(value) && any(vapply(value, contains_environment, logical(1))))
  }
  run_named_test("cached fits report complete loading, preparation and fitting stages", {
    cache <- make_cache("whole")
    events <- character()
    observer <- function(stage, event) events <<- c(events, paste(stage, event, sep = ":"))
    tryCatch({
      finemapPgenPhenotype("left", se, reader, covariates, distance, "NA",
        timing_callback = observer, genotype_cache = cache)
      expect_identical_value(events, c("data_loading:start", "data_loading:complete",
        "data_preparation:start", "data_preparation:complete",
        "fine_mapping:start", "fine_mapping:complete"), "cached stage events")
    }, finally = clearPgenGenotypeCache(cache))
  })
  run_named_test("explicit no-cache input preserves the current native-first fit", {
    reference <- fit("left", explicit = FALSE)
    compare_fit(fit("left"), reference, exact = TRUE)
  })

  references <- stats::setNames(lapply(phenotype_ids, fit), phenotype_ids)
  run_named_test("cache fit comparison separates optimizer rounding from scientific changes", {
    reference <- references$left
    rounded <- reference
    mu_index <- which.max(abs(reference$result$mu))
    pip_index <- which.max(reference$result$pip)
    rounded$result$mu[[mu_index]] <- rounded$result$mu[[mu_index]] + 1e-8
    rounded$result$pip[[pip_index]] <- rounded$result$pip[[pip_index]] - 2e-9
    compare_fit(rounded, reference)
    changed <- reference
    changed$result$pip[[pip_index]] <- changed$result$pip[[pip_index]] - 1e-5
    expect_error_message(compare_fit(changed, reference), "PIP absolute")
    changed <- reference
    changed$result$sets$cs <- list(L1 = 999L)
    expect_error_message(compare_fit(changed, reference), "credible-set memberships")
    changed <- reference
    changed$result$converged <- !changed$result$converged
    expect_error_message(compare_fit(changed, reference), "convergence status")
    changed <- reference
    changed$result$niter <- changed$result$niter + 1L
    expect_error_message(compare_fit(changed, reference), "iteration count")
    changed <- reference
    changed$result$variant_id <- rev(changed$result$variant_id)
    expect_error_message(compare_fit(changed, reference), "fit variant order")
    changed <- reference
    changed$rng[[1L]] <- changed$rng[[1L]] + 1L
    expect_error_message(compare_fit(changed, reference), "fit RNG state")
    changed <- reference
    changed$inputs$X[[1L]] <- changed$inputs$X[[1L]] + 1e-9
    expect_error_message(compare_fit(changed, reference), "Adjusted fit input X")
    expect_error_message(compare_fit(rounded, reference, exact = TRUE), "guarded complete fit")
  })
  run_named_test("whole and rolling caches preserve overlapping, disjoint and returning fits", {
    regions <- tibble::tibble(phenotype_id = phenotype_ids, region = worklist$region)
    reference_fits <- lapply(references, `[[`, "result")
    reference_tables <- formatSusieOutputs(reference_fits, regions)
    for (mode in c("whole", "rolling")) {
      cache <- make_cache(mode)
      tryCatch({
        actual <- stats::setNames(lapply(phenotype_ids, function(id) fit(id, cache = cache)),
                                  phenotype_ids)
        for (id in phenotype_ids) compare_fit(actual[[id]], references[[id]])
        actual_fits <- lapply(actual, `[[`, "result")
        tables <- formatSusieOutputs(actual_fits, regions)
        comparison <- all.equal(tables, reference_tables, tolerance = model_tolerance,
                                check.attributes = TRUE)
        if (!isTRUE(comparison)) stop(paste(comparison, collapse = "\n"), call. = FALSE)
        saved <- file.path(directory, paste0(mode, ".rds"))
        saveRDS(actual_fits, saved)
        stopifnot(!contains_environment(readRDS(saved)))
        expect_identical_value(lapply(actual_fits, `[[`, "variant_id"),
                               lapply(reference_fits, `[[`, "variant_id"))
      }, finally = clearPgenGenotypeCache(cache))
    }
  })
  run_named_test("prepared fits reject changed samples, covariates and feature windows", {
    cache <- make_cache("whole")
    tryCatch({
      changed_covariates <- covariates
      changed_covariates[selected_samples[[1L]], "pc1"] <-
        changed_covariates[selected_samples[[1L]], "pc1"] + 0.1
      expect_error_message(fit("left", cov = changed_covariates, cache = cache),
                           "sample order or covariates differ")
      reordered <- se[, rev(seq_len(ncol(se)))]
      expect_error_message(finemapPgenPhenotype("left", reordered, reader, covariates,
        distance, "NA", genotype_cache = cache), "sample order or covariates differ")
      expect_error_message(finemapPgenPhenotype("left", se, reader, covariates,
        distance, "NA", genotype_cache = cache, cache_outcome_id = "distant"),
        "window differs")
      stopifnot(pgenGenotypeCacheSummary(cache)$variants_read == 0)
    }, finally = clearPgenGenotypeCache(cache))
  })
  run_named_test("genotype-triggered dense guards preserve exact complete fits", {
    guarded_cov <- cbind(genotype = genotypes[7L, selected_samples],
                          pc = covariates[selected_samples, "pc2"])
    rownames(guarded_cov) <- selected_samples
    reference <- fit("left", cov = guarded_cov)
    stopifnot(variant_ids[[7L]] %in%
                reference$result$covariate_explained_variants$variant_id)
    for (mode in c("whole", "rolling")) {
      cache <- make_cache(mode, guarded_cov)
      events <- character()
      observer <- function(stage, event) events <<- c(events, paste(stage, event, sep = ":"))
      tryCatch({
        compare_fit(fit("left", cov = guarded_cov, cache = cache,
                        timing_callback = observer), reference, exact = TRUE)
        stopifnot(!any(grepl(":failed$", events)))
        for (stage in c("data_loading", "data_preparation", "fine_mapping")) {
          stopifnot(sum(events == paste0(stage, ":start")) ==
                      sum(events == paste0(stage, ":complete")))
        }
        stopifnot(sum(events == "data_loading:start") == 2L,
                  sum(events == "data_preparation:start") == 2L,
                  sum(events == "fine_mapping:start") == 1L)
      }, finally = clearPgenGenotypeCache(cache))
    }
  })
  run_named_test("phenotype-triggered dense guards preserve exact complete fits", {
    y <- phenotype_values["left", ]
    normalized <- qnorm((rank(y) - 0.5) / length(y))
    guarded_cov <- cbind(response = normalized + 1e-10 * sin(seq_along(y)),
                          pc = covariates[selected_samples, "pc2"])
    rownames(guarded_cov) <- selected_samples
    reference <- fit("left", cov = guarded_cov)
    for (mode in c("whole", "rolling")) {
      cache <- make_cache(mode, guarded_cov)
      tryCatch(compare_fit(fit("left", cov = guarded_cov, cache = cache), reference,
                           exact = TRUE), finally = clearPgenGenotypeCache(cache))
    }
  })
  run_named_test("cached multi-outcome gene artifacts contain fits but no cache state", {
    outcomes <- stats::setNames(lapply(phenotype_ids[1:2], function(id) {
      list(phenotype_id = id, se = se, covariates = covariates,
           cis_distance = distance, region = worklist$region[match(id, phenotype_ids)])
    }), phenotype_ids[1:2])
    reference <- runPgenGene("paired", outcomes, reader,
      file.path(directory, "baseline"), 341L, "NA")
    cache <- make_cache("rolling")
    tryCatch({
      actual <- runPgenGene("paired", outcomes, reader,
        file.path(directory, "cached"), 341L, "NA", genotype_cache = cache)
      timing_columns <- names(actual$index)[grepl("_seconds$", names(actual$index))]
      required_timing <- c("data_loading_seconds", "data_preparation_seconds",
                           "fine_mapping_seconds", "outcome_total_seconds")
      stopifnot(all(required_timing %in% timing_columns),
                identical(names(actual$index), names(reference$index)))
      scientific_columns <- setdiff(names(actual$index), timing_columns)
      expect_identical_value(actual$index[scientific_columns],
                             reference$index[scientific_columns], "scientific multi-outcome index")
      for (index in list(actual$index, reference$index)) {
        timings <- index[timing_columns]
        stopifnot(all(vapply(timings, is.numeric, logical(1))),
                  all(is.finite(as.matrix(timings))), all(as.matrix(timings) >= 0),
                  all(index$outcome_total_seconds + 1e-9 >=
                    rowSums(as.matrix(index[setdiff(timing_columns, "outcome_total_seconds")]))))
      }
      for (name in c("summary", "credible_sets", "diagnostics")) {
        comparison <- all.equal(actual[[name]], reference[[name]], tolerance = model_tolerance,
                                check.attributes = TRUE)
        if (!isTRUE(comparison)) stop(paste(comparison, collapse = "\n"), call. = FALSE)
      }
      saved <- readRDS(actual$files[["rds"]])
      stopifnot(!contains_environment(saved))
      expected <- readRDS(reference$files[["rds"]])
      expect_identical_value(names(saved), names(expected), "saved outcome order")
      for (id in names(expected)) compare_model(saved[[id]], expected[[id]])
    }, finally = clearPgenGenotypeCache(cache))
  })
  message("All prepared PGEN cache fit tests passed.")
  invisible(TRUE)
}
