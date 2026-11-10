# Compare the sufficient-statistics path with real individual-data fits.
run_susie_sufficient_statistics_tests <- function(repo_dir = ".") {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  source(file.path(repo_dir, "tests", "test_helpers.R"), local = TRUE)
  functions_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH",
                               unset = file.path(repo_dir, "R", "utils"))
  for (filename in c("SusiePerformance.R", "SusieFunctions.R")) {
    source(file.path(functions_path, filename), local = TRUE)
  }
  set.seed(7219L)
  n <- 200L
  p <- 20L
  X <- matrix(rnorm(n * p), nrow = n,
    dimnames = list(paste0("s", seq_len(n)), paste0("v", seq_len(p))))
  # Offsets and unequal scales expose a missing centering or scale correction.
  X <- sweep(sweep(X, 2L, seq(0.5, 2, length.out = p), "*"),
             2L, seq(-2, 2, length.out = p), "+")
  responses <- list(
    strong = 1.3 * X[, 2L] - 0.8 * X[, 14L] + rnorm(n, sd = .6),
    weak = .2 * X[, 8L] + rnorm(n)
  )
  statistics <- susieR::compute_suff_stat(X, responses$strong)[
    c("XtX", "X_colmeans", "n")]
  capture_fit <- function(y, sufficient_statistics = NULL, legacy = FALSE) {
    set.seed(1841L)
    output <- capture.output(value <- suppressMessages(fitSusieModel(
      X, cbind(y), legacy = legacy,
      sufficient_statistics = sufficient_statistics
    )))
    list(fit = value, rng = .Random.seed)
  }
  assert_close <- function(actual, expected, label, tolerance = 1e-8) {
    if (!identical(attributes(actual), attributes(expected)) ||
        !is.numeric(actual) || !is.numeric(expected) ||
        length(actual) != length(expected) ||
        any(!is.finite(actual)) || any(!is.finite(expected))) {
      stop(label, " has a different shape or nonfinite values", call. = FALSE)
    }
    difference <- max(c(0, abs(actual - expected)))
    if (difference > tolerance) {
      stop(label, " maximum absolute difference is ",
           sprintf("%.17g", difference), call. = FALSE)
    }
    invisible(TRUE)
  }
  compare_science <- function(actual, expected) {
    assert_close(actual$fit$pip, expected$fit$pip, "PIP")
    assert_close(susieR::susie_get_posterior_mean(actual$fit),
                 susieR::susie_get_posterior_mean(expected$fit), "posterior mean")
    assert_close(susieR::susie_get_posterior_sd(actual$fit),
                 susieR::susie_get_posterior_sd(expected$fit), "posterior SD")
    assert_close(actual$fit$z, expected$fit$z, "marginal z")
    assert_close(actual$fit$Xr, expected$fit$Xr, "raw-scale fitted component")
    assert_close(actual$fit$fitted, expected$fit$fitted, "fitted response")
    assert_close(actual$fit$intercept, expected$fit$intercept, "intercept")
    expect_identical_value(actual$fit$sets$cs, expected$fit$sets$cs, "CS membership")
    expect_identical_value(actual$fit$sets$cs_index, expected$fit$sets$cs_index,
                           "CS effects")
    if (!isTRUE(all.equal(actual$fit$sets$purity, expected$fit$sets$purity,
                         tolerance = 1e-12))) {
      stop("CS purity differs from native individual-data purity", call. = FALSE)
    }
    expect_identical_value(actual$fit$niter, expected$fit$niter, "fit iterations")
    expect_identical_value(actual$fit$converged, expected$fit$converged, "fit convergence")
    expect_identical_value(actual$rng, expected$rng, "fit RNG state")
    expect_identical_value(actual$fit$covariate_explained_variants,
      expected$fit$covariate_explained_variants, "covariate diagnostics")
  }
  run_named_test("sufficient statistics preserve native scientific fit values", {
    # This direct call also fails when the new argument has not been added.
    first <- capture_fit(responses$strong, statistics)
    stopifnot("The pinned source must enable sufficient-statistics fits" =
      isTRUE(.susieSufficientStatsSourceSupported(asNamespace("susieR"))))
    reference <- capture_fit(responses$strong)
    stopifnot("The strong fixture must produce a credible set" =
      length(reference$fit$sets$cs) > 0L)
    compare_science(first, reference)
    weak_reference <- capture_fit(responses$weak)
    observer <- new.env(parent = emptyenv())
    observer$calls <- 0L
    tracer <- substitute({
      .ss_fit_observer <- captured
      .ss_fit_observer$calls <- .ss_fit_observer$calls + 1L
    }, list(captured = observer))
    suppressMessages(trace("susie_ss", tracer = tracer,
      where = asNamespace("susieR"), print = FALSE))
    # Tracing changes the audited body. Permit this one real, traced call only.
    saved_guard <- .susieSufficientStatsSourceSupported
    .susieSufficientStatsSourceSupported <- function(ns) TRUE
    tryCatch({
      compare_science(capture_fit(responses$weak, statistics), weak_reference)
      stopifnot("Eligible statistics must use the real native SS fitter" =
        observer$calls == 1L)
    }, finally = {
      .susieSufficientStatsSourceSupported <- saved_guard
      suppressMessages(untrace("susie_ss", where = asNamespace("susieR")))
    })
  })
  run_named_test("legacy fits ignore sufficient statistics", {
    expect_identical_value(capture_fit(responses$strong, statistics, legacy = TRUE),
      capture_fit(responses$strong, legacy = TRUE), "legacy complete native fit")
  })
  run_named_test("summary fit finalization keeps native CS purity defaults", {
    reference <- capture_fit(responses$strong)
    observer <- new.env(parent = emptyenv())
    observer$calls <- 0L
    tracer <- substitute({
      .ss_cs_observer <- captured
      .ss_cs_observer$calls <- .ss_cs_observer$calls + 1L
      .ss_cs_observer$X <- X
      .ss_cs_observer$Xcorr <- Xcorr
      .ss_cs_observer$n_purity <- n_purity
    }, list(captured = observer))
    suppressMessages(trace("susie_get_cs", tracer = tracer,
      where = asNamespace("susieR"), print = FALSE))
    tryCatch({
      pending <- reference$fit
      pending$sets <- NULL
      set.seed(1841L)
      seed_before <- .Random.seed
      finalized <- finalizeSusieSufficientStatistics(pending, X,
        cbind(responses$strong), statistics)
      stopifnot("Finalization must calculate CS once with native defaults" =
        observer$calls == 1L)
      expect_identical_value(observer$X, X, "CS individual-data input")
      expect_identical_value(observer$Xcorr, NULL, "CS correlation input")
      expect_identical_value(observer$n_purity, "auto", "CS purity default")
      expect_identical_value(finalized$sets, reference$fit$sets, "native CS result")
      expect_identical_value(.Random.seed, seed_before, "CS finalization RNG state")
    }, finally = suppressMessages(untrace("susie_get_cs",
      where = asNamespace("susieR"))))
  })
  run_named_test("summary fit finalization preserves native purity sampling RNG", {
    broad_X <- X[, rep(1L, p)] +
      matrix(1e-3 * sin(seq_len(n * p)), nrow = n)
    dimnames(broad_X) <- dimnames(X)
    broad_statistics <- susieR::compute_suff_stat(broad_X, responses$strong)[
      c("XtX", "X_colmeans", "n")]
    pending <- reference$fit
    pending$alpha[,] <- 1 / p
    pending$V[] <- 0
    pending$V[[1L]] <- 1
    # Force a small cap only in this test. Native sample() still selects the
    # purity subset, so a small matrix can test actual RNG use.
    suppressMessages(trace("get_purity", tracer = quote(n <- 2L),
      where = asNamespace("susieR"), print = FALSE))
    tryCatch({
      set.seed(2887L)
      seed_before <- .Random.seed
      expected <- susieR::susie_get_cs(pending, X = broad_X,
        coverage = .95, min_abs_corr = .5, n_purity = "auto")
      expected_rng <- .Random.seed
      stopifnot("The broad fixture must retain a CS and sample its purity" =
        length(expected$cs) > 0L && !identical(expected_rng, seed_before))
      set.seed(2887L)
      finalized <- finalizeSusieSufficientStatistics(pending, broad_X,
        cbind(responses$strong), broad_statistics)
      expect_identical_value(finalized$sets, expected, "sampled native CS result")
      expect_identical_value(.Random.seed, expected_rng, "sampled native CS RNG state")
    }, finally = suppressMessages(untrace("get_purity",
      where = asNamespace("susieR"))))
  })
  run_named_test("changed native SS calculations use the exact individual-data fit", {
    namespace <- asNamespace("susieR")
    stopifnot("The unchanged native SS source must pass its guard" =
      isTRUE(.susieSufficientStatsSourceSupported(namespace)))
    reference <- capture_fit(responses$strong)
    original <- get("get_ER2.ss", envir = namespace, inherits = FALSE)
    was_locked <- bindingIsLocked("get_ER2.ss", namespace)
    tryCatch({
      if (was_locked) unlockBinding("get_ER2.ss", namespace)
      assign("get_ER2.ss", function(data, model) 0, envir = namespace)
      lockBinding("get_ER2.ss", namespace)
      stopifnot("A changed native SS calculation must fail its source guard" =
        identical(.susieSufficientStatsSourceSupported(namespace), FALSE))
      stopifnot("An SS-only change must preserve the individual-data source guard" =
        isTRUE(.susieFitSourceSupported(namespace)))
      expect_identical_value(capture_fit(responses$strong, statistics), reference,
                             "changed-source complete native fallback")
    }, finally = {
      if (bindingIsLocked("get_ER2.ss", namespace)) {
        unlockBinding("get_ER2.ss", namespace)
      }
      assign("get_ER2.ss", original, envir = namespace)
      if (was_locked) lockBinding("get_ER2.ss", namespace)
    })
  })
  run_named_test("unsupported sufficient-statistics source uses the exact native fit", {
    reference <- capture_fit(responses$strong)
    saved_guard <- .susieSufficientStatsSourceSupported
    .susieSufficientStatsSourceSupported <- function(ns) FALSE
    tryCatch({
      expect_identical_value(capture_fit(responses$strong, statistics), reference,
                             "unsupported-source complete native fit")
    }, finally = .susieSufficientStatsSourceSupported <- saved_guard)
  })
  run_named_test("sufficient statistics reject a different sample count", {
    invalid <- statistics
    invalid$n <- n - 1L
    expect_error_message(capture_fit(responses$strong, invalid), "sample")
  })
  run_named_test("sufficient statistics reject a different matrix dimension", {
    invalid <- statistics
    invalid$XtX <- invalid$XtX[-1L, -1L, drop = FALSE]
    expect_error_message(capture_fit(responses$strong, invalid), "dimension")
  })
  run_named_test("sufficient statistics reject a different variant order", {
    invalid <- statistics
    dimnames(invalid$XtX) <- rep(list(rev(colnames(X))), 2L)
    expect_error_message(capture_fit(responses$strong, invalid), "variant")
  })
  message("All sufficient-statistics fit tests passed.")
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  benchmark_library <- Sys.getenv("BENCH_LIBRARY", unset = "")
  if (nzchar(benchmark_library)) .libPaths(c(benchmark_library, .libPaths()))
  args <- commandArgs(trailingOnly = TRUE)
  repo_dir <- normalizePath(if (length(args)) args[[1L]] else getwd(), mustWork = TRUE)
  run_susie_sufficient_statistics_tests(repo_dir)
}
