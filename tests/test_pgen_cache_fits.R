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

  # Keep the exact inputs at the real fitter boundary. Do not change the fit.
  original_fit_susie_model <- fitSusieModel
  fit_inputs <- NULL
  fitSusieModel <- function(X, y, ...) {
    fit_inputs <<- list(X = X, y = y)
    original_fit_susie_model(X, y, ...)
  }
  fit <- function(id, cov = covariates, cache = NULL, explicit = TRUE,
                   timing_callback = NULL) {
    set.seed(259L)
    fit_inputs <<- NULL
    arguments <- list(phenotype_id = id, se = se, reader = reader, covariates = cov,
                      cis_distance = distance, missing_code = "NA")
    if (explicit) {
      arguments["genotype_cache"] <- list(cache)
      arguments$cache_outcome_id <- id
    }
    if (!is.null(timing_callback)) arguments$timing_callback <- timing_callback
    result <- do.call(finemapPgenPhenotype, arguments)
    list(result = result, rng = .Random.seed, inputs = fit_inputs, phenotype_id = id)
  }
  make_cache <- function(mode, cov = covariates) {
    plan <- planPgenGenotypeCache(reader, worklist, selected_samples, mode = mode,
      cache_memory_gb = 0.1, job_memory_gb = 1)
    createPgenGenotypeCache(reader, plan, selected_samples,
      cbind(cov[selected_samples, , drop = FALSE], 1))
  }
  report_fit_differences <- function(actual, reference) {
    max_absolute_error <- function(actual, reference) {
      if (!is.numeric(actual) || !is.numeric(reference) ||
          length(actual) != length(reference) ||
          !identical(dim(actual), dim(reference))) return(NA_real_)
      if (!length(actual)) return(0)
      max(abs(actual - reference))
    }
    dimensions <- function(value) {
      if (is.null(dim(value))) "vector" else paste(dim(value), collapse = "x")
    }
    report_field <- function(field, actual, reference) {
      cat("CACHE FIT DIAGNOSTIC field=", field,
          " actual_length=", length(actual),
          " reference_length=", length(reference),
          " actual_dimensions=", dimensions(actual),
          " reference_dimensions=", dimensions(reference),
          " max_abs_error=", sprintf("%.17g", max_absolute_error(actual, reference)),
          " actual_max_abs=", sprintf("%.17g", max(abs(actual))),
          " reference_max_abs=", sprintf("%.17g", max(abs(reference))),
          "\n", sep = "")
    }
    cat("CACHE FIT DIAGNOSTIC phenotype=", actual$phenotype_id,
        " actual_converged=", actual$result$converged,
        " reference_converged=", reference$result$converged,
        " credible_sets_identical=",
        identical(actual$result$sets$cs, reference$result$sets$cs), "\n", sep = "")
    for (field in c("X", "y")) {
      report_field(paste0("input_", field), actual$inputs[[field]], reference$inputs[[field]])
    }
    for (field in c("alpha", "mu", "mu2", "V", "lbf", "lbf_variable", "pip")) {
      report_field(field, actual$result[[field]], reference$result[[field]])
    }
  }
  # Linux OpenBLAS produced a 2.22e-16 input difference with identical y,
  # convergence and credible sets. The native Brent prior optimizer uses
  # R's default sqrt(machine epsilon) tolerance, about 1.49e-8. It amplified
  # the input error to a 6.84e-10 PIP and 6.79e-8 log-BF difference. Bound inputs and
  # native scientific values separately from internal optimizer roundoff.
  # Guarded fits, identities, RNG state and convergence must remain exact.
  input_absolute_tolerance <- 1e-13
  scientific_absolute_tolerance <- 1e-8
  fit_tolerance <- 1e-7
  compare_absolute <- function(actual, reference, tolerance, label) {
    same_structure <- is.numeric(actual) && is.numeric(reference) &&
      length(actual) == length(reference) &&
      identical(attributes(actual), attributes(reference)) &&
      all(is.finite(actual)) && all(is.finite(reference))
    error <- if (same_structure) max(c(0, abs(actual - reference))) else Inf
    if (!is.finite(error) || error > tolerance) {
      stop(label, " max absolute error exceeds ", sprintf("%.1e", tolerance),
           ": ", sprintf("%.17g", error), call. = FALSE)
    }
    invisible(TRUE)
  }
  compare_model_values <- function(actual, reference, exact = FALSE) {
    expect_identical_value(actual$variant_id, reference$variant_id, "fit variant order")
    expect_identical_value(actual$sets$cs, reference$sets$cs, "credible-set memberships")
    expect_identical_value(actual$sets$cs_index, reference$sets$cs_index, "credible-set effects")
    expect_identical_value(actual$converged, reference$converged, "fit convergence status")
    if (exact) {
      expect_identical_value(actual, reference, "guarded complete fit")
      return(invisible(TRUE))
    }
    compare_absolute(actual$pip, reference$pip,
                     scientific_absolute_tolerance, "native PIP")
    compare_absolute(susieR::susie_get_posterior_mean(actual),
                     susieR::susie_get_posterior_mean(reference),
                     scientific_absolute_tolerance, "native posterior mean")
    compare_absolute(susieR::susie_get_posterior_sd(actual),
                     susieR::susie_get_posterior_sd(reference),
                     scientific_absolute_tolerance, "native posterior SD")
    compare_absolute(actual$z, reference$z, scientific_absolute_tolerance, "marginal z")
    # Prior variance is an optimized hyperparameter. Validate it with the same
    # explicit component bound as its exported table, then avoid a second,
    # relative-only check that is tighter on different BLAS implementations.
    compare_output_values(tibble::tibble(V = actual$V), tibble::tibble(V = reference$V))
    comparison_actual <- actual
    comparison_actual$V <- reference$V
    comparison <- all.equal(comparison_actual, reference, tolerance = fit_tolerance,
                            check.attributes = TRUE)
    if (!isTRUE(comparison)) stop(paste(comparison, collapse = "\n"), call. = FALSE)
    invisible(TRUE)
  }
  compare_fit <- function(actual, reference, exact = FALSE) {
    tryCatch({
      for (field in c("X", "y")) {
        compare_absolute(actual$inputs[[field]], reference$inputs[[field]],
                         input_absolute_tolerance, paste("fitter input", field))
      }
      expect_identical_value(actual$rng, reference$rng, "fit RNG state")
      compare_model_values(actual$result, reference$result, exact = exact)
    }, error = function(error) {
      report_fit_differences(actual, reference)
      stop(conditionMessage(error), call. = FALSE)
    })
  }
  compare_output_values <- function(actual, reference) {
    actual_tables <- if (is.data.frame(actual)) list(actual) else actual
    reference_tables <- if (is.data.frame(reference)) list(reference) else reference
    expect_identical_value(names(actual), names(reference), "output names")
    expect_identical_value(attributes(actual), attributes(reference), "output attributes")
    expect_identical_value(length(actual_tables), length(reference_tables), "output table count")
    key_columns <- c("molecular_trait_id", "phenotype_id", "outcome_id", "gene_id",
      "variant", "variant_id", "chromosome", "position", "pos", "ref", "alt",
      "cs_id", "cs_index", "region", "cs_size", "sample_id", "sample_ids",
      "sample_index", "sample_indices", "variant_index", "variant_indices", "sample_count")
    for (i in seq_along(reference_tables)) {
      actual_table <- actual_tables[[i]]
      reference_table <- reference_tables[[i]]
      table_name <- if (is.null(names(reference_tables))) as.character(i) else
        names(reference_tables)[[i]]
      label <- paste("output", table_name)
      expect_identical_value(attributes(actual_table), attributes(reference_table),
                             paste(label, "table attributes"))
      scientific_columns <- intersect(c("pip", "posterior_mean", "posterior_sd", "z",
                                         "X_column_scale_factors"),
                                      names(reference_table))
      for (field in scientific_columns) {
        compare_absolute(actual_table[[field]], reference_table[[field]],
                         scientific_absolute_tolerance, paste(label, field))
      }
      for (field in names(reference_table)) {
        a <- actual_table[[field]]
        b <- reference_table[[field]]
        column_label <- paste(label, field)
        expect_identical_value(typeof(a), typeof(b), paste(column_label, "type"))
        expect_identical_value(attributes(a), attributes(b),
                               paste(column_label, "attributes"))
        if (!is.numeric(b) || field %in% key_columns ||
            grepl("(^|_)(id|ids|index|indices)$", field)) {
          expect_identical_value(a, b, column_label)
          next
        }
        expect_identical_value(is.na(a), is.na(b), paste(column_label, "NA mask"))
        expect_identical_value(is.nan(a), is.nan(b), paste(column_label, "NaN mask"))
        expect_identical_value(is.infinite(a), is.infinite(b),
                               paste(column_label, "infinite mask"))
        expect_identical_value(a[is.infinite(a)], b[is.infinite(b)],
                               paste(column_label, "infinite values"))
        finite <- is.finite(b)
        difference <- abs(a[finite] - b[finite])
        # Scientific fields already passed their stricter absolute check above.
        # Use the existing 1e-7 fit tolerance as both an absolute and relative
        # bound for components: near-zero log-BFs differ by about 5e-8 on Linux.
        allowed <- fit_tolerance + fit_tolerance * abs(b[finite])
        if (any(difference > allowed)) {
          stop(column_label, " max absolute error ",
               sprintf("%.17g", max(difference)),
               " exceeds the absolute and relative bound", call. = FALSE)
        }
      }
    }
    invisible(TRUE)
  }
  contains_environment <- function(value) {
    is.environment(value) ||
      (is.list(value) && any(vapply(value, contains_environment, logical(1))))
  }
  run_named_test("output comparison accepts negligible near-zero component roundoff", {
    reference <- list(lbf_df = tibble::tibble(lbf_variable3 = 1e-4))
    actual <- list(lbf_df = tibble::tibble(lbf_variable3 = 1e-4 + 5.382028e-8))
    compare_output_values(actual, reference)
  })
  run_named_test("output comparison rejects material near-zero component changes", {
    for (field in c("lbf_variable3", "mu_3")) {
      for (change in c(2e-7, 1e-6)) {
        reference <- list(component = tibble::as_tibble(stats::setNames(list(0), field)))
        actual <- list(component = tibble::as_tibble(stats::setNames(list(change), field)))
        expect_error_message(compare_output_values(actual, reference), field)
      }
    }
  })
  run_named_test("output comparison preserves scientific absolute bounds", {
    for (field in c("pip", "posterior_mean", "posterior_sd", "z", "X_column_scale_factors")) {
      reference <- list(science = tibble::as_tibble(stats::setNames(list(.9), field)))
      actual <- list(science = tibble::as_tibble(stats::setNames(list(.9 + 5e-8), field)))
      expect_error_message(compare_output_values(actual, reference), field)
    }
  })
  run_named_test("output comparison requires exact output keys", {
    reference <- list(variant_df = tibble::tibble(variant = "v1", position = 1e8,
                                                 cs_index = 1))
    for (field in names(reference$variant_df)) {
      actual <- reference
      actual$variant_df[[field]][[1L]] <- if (field == "variant") "v2" else
        reference$variant_df[[field]][[1L]] + 1
      expect_error_message(compare_output_values(actual, reference), field)
    }
  })
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
  run_named_test("cache model comparison bounds prior optimizer roundoff", {
    reference <- references$left$result
    rounded <- reference
    rounded$V[[1L]] <- rounded$V[[1L]] + 1.331136e-7
    compare_model_values(rounded, reference)
    changed <- reference
    changed$V[[1L]] <- changed$V[[1L]] + 1e-6
    expect_error_message(compare_model_values(changed, reference), "V")
  })
  run_named_test("overlapping whole-cache windows always use the ordinary fitter", {
    nearby <- worklist[1:2, , drop = FALSE]
    plan <- planPgenGenotypeCache(reader, nearby, selected_samples, mode = "whole",
      cache_memory_gb = 0.1, job_memory_gb = 1)
    cache <- createPgenGenotypeCache(reader, plan, selected_samples,
      cbind(covariates[selected_samples, , drop = FALSE], 1))
    original_ordinary_fitter <- runSusieWithFitOptimizations
    ordinary_fit_calls <- 0L
    runSusieWithFitOptimizations <- function(...) {
      ordinary_fit_calls <<- ordinary_fit_calls + 1L
      original_ordinary_fitter(...)
    }
    tryCatch({
      for (id in phenotype_ids[1:2]) {
        previous_calls <- ordinary_fit_calls
        actual <- fit(id, cache = cache)
        expect_identical_value(ordinary_fit_calls, previous_calls + 1L,
                               "ordinary fitter calls for each cached outcome")
        compare_fit(actual, references[[id]])
      }
      stopifnot(!"sufficient_statistics" %in% ls(cache),
                !any(vapply(mget(ls(cache), envir = cache, inherits = FALSE), function(value)
                  is.list(value) && "XtX" %in% names(value), logical(1))))
    }, finally = {
      runSusieWithFitOptimizations <- original_ordinary_fitter
      clearPgenGenotypeCache(cache)
    })
  })
  run_named_test("cache fit comparison rejects changed fitter inputs", {
    reference <- references[["overlap"]]
    for (field in c("X", "y")) {
      changed <- reference
      changed$inputs[[field]][[1L]] <- changed$inputs[[field]][[1L]] + 1e-6
      invisible(capture.output(expect_error_message(
        compare_fit(changed, reference), paste("fitter input", field)
      )))
    }
  })
  run_named_test("cache fit comparison rejects a material native PIP change", {
    reference <- references[["overlap"]]
    changed <- reference
    changed$result$pip[[1L]] <- changed$result$pip[[1L]] + 1e-6
    invisible(capture.output(expect_error_message(
      compare_fit(changed, reference), "native PIP"
    )))
  })
  run_named_test("cache fit comparison accepts machine-precision input changes", {
    reference <- references[["overlap"]]
    rounded <- reference
    for (field in c("X", "y")) {
      rounded$inputs[[field]][[1L]] <-
        rounded$inputs[[field]][[1L]] + .Machine$double.eps
    }
    compare_fit(rounded, reference)
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
        compare_output_values(tables, reference_tables)
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
        compare_output_values(actual[[name]], reference[[name]])
      }
      saved <- readRDS(actual$files[["rds"]])
      stopifnot(!contains_environment(saved))
      expected <- readRDS(reference$files[["rds"]])
      expect_identical_value(names(saved), names(expected), "saved outcome order")
      for (id in names(expected)) compare_model_values(saved[[id]], expected[[id]])
    }, finally = clearPgenGenotypeCache(cache))
  })
  message("All prepared PGEN cache fit tests passed.")
  invisible(TRUE)
}
