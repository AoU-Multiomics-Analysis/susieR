# Real fits check that transient input statistics preserve all model outputs.
run_marginal_stat_reuse_tests <- function(repo_dir = ".") {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  strict_same <- function(x, y) {
    identical(x, y, num.eq = FALSE, single.NA = FALSE, attrib.as.set = FALSE)
  }
  expect_same <- function(expected, actual, label) {
    if (!strict_same(expected, actual)) stop(label, call. = FALSE)
    invisible(TRUE)
  }
  functions_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH",
                              unset = file.path(repo_dir, "R", "utils"))
  make_helpers <- function() {
    helpers <- new.env(parent = globalenv())
    sys.source(file.path(functions_path, "SusiePerformance.R"), envir = helpers)
    sys.source(file.path(functions_path, "SusieFunctions.R"), envir = helpers)
    helpers
  }
  helpers <- make_helpers()
  native <- getFromNamespace("susie", "susieR")
  namespace <- asNamespace("susieR")
  namespace_names <- ls(namespace, all.names = TRUE)
  namespace_values <- mget(namespace_names, envir = namespace, inherits = FALSE)
  namespace_locks <- vapply(namespace_names, bindingIsLocked, logical(1), env = namespace)
  check_namespace <- function() {
    expect_same(namespace_names, ls(namespace, all.names = TRUE), "Native names must not change")
    expect_same(namespace_values, mget(namespace_names, envir = namespace, inherits = FALSE),
                "Native functions and values must not change")
    expect_same(namespace_locks,
                vapply(namespace_names, bindingIsLocked, logical(1), env = namespace),
                "Native binding locks must not change")
  }

  set.seed(76321)
  X <- matrix(rnorm(48L * 32L), nrow = 48L,
              dimnames = list(paste0("s", 1:48), paste0("v", 1:32)))
  y <- .5 * X[, 1L] - .3 * X[, 9L] + rnorm(48L)
  callback_calls <- 0L
  # Removing callback support must fail on a real native fit.
  callback_fit <- helpers$runSusieWithFitOptimizations(X, y, L = 3L, max_iter = 10L,
    .postprocess = function(fit, stats) {
      callback_calls <<- callback_calls + 1L
      stopifnot("A callback must receive a completed native fit" = inherits(fit, "susie"))
      fit
    })
  stopifnot("Post-fit processing must run exactly once" = callback_calls == 1L,
            inherits(callback_fit, "susie"))

  # These spies still call the real calculations. They count matrix scans
  # without tracing or changing a package namespace.
  counts <- new.env(parent = emptyenv())
  reset_counts <- function() {
    for (name in c("mean", "sd", "range", "native_z", "fit_mean", "fit_sd")) {
      counts[[name]] <- 0L
    }
  }
  reset_counts()
  tracked_X <- X
  helpers$colMeans <- function(x, ...) {
    if (strict_same(x, tracked_X)) counts$mean <- counts$mean + 1L
    base::colMeans(x, ...)
  }
  helpers$range <- function(..., na.rm = FALSE) {
    inputs <- list(...)
    if (length(inputs) == 1L && strict_same(inputs[[1L]], tracked_X)) {
      counts$range <- counts$range + 1L
    }
    base::range(..., na.rm = na.rm)
  }
  helpers$`::` <- function(package, name) {
    package_name <- as.character(substitute(package))
    function_name <- as.character(substitute(name))
    calculation <- getExportedValue(package_name, function_name)
    if (identical(package_name, "matrixStats") && identical(function_name, "colSds")) {
      return(function(x, ...) {
        if (strict_same(x, tracked_X)) counts$sd <- counts$sd + 1L
        calculation(x, ...)
      })
    }
    calculation
  }
  helpers$getFromNamespace <- function(name, package) {
    calculation <- utils::getFromNamespace(name, package)
    if (identical(name, "calc_z") && identical(package, "susieR")) {
      return(function(X, Y, center = FALSE, scale = FALSE) {
        counts$native_z <- counts$native_z + 1L
        calculation(X, Y, center = center, scale = scale)
      })
    }
    calculation
  }
  make_stats <- function(x) {
    list(X = x, n = nrow(x), p = ncol(x), finite = TRUE,
         x_mean = base::colMeans(x), x_sd = matrixStats::colSds(x))
  }
  stats <- make_stats(X)
  ordinary_z <- helpers$computeMarginalZScores(X, y)
  stopifnot("The control must exercise each old matrix scan" =
              counts$mean == 1L && counts$sd == 1L && counts$range == 1L,
            counts$native_z == 0L)
  reset_counts()
  expect_same(ordinary_z, helpers$computeMarginalZScores(X, y, .x_stats = stats),
              "Reused statistics must preserve every marginal-score bit")
  stopifnot("Valid statistics must remove all three repeated matrix scans" =
              counts$mean == 0L && counts$sd == 0L && counts$range == 0L,
            counts$native_z == 0L)
  for (response in list(y, cbind(y), cbind(first = y, second = y + X[, 5L]))) {
    for (center in c(FALSE, TRUE)) {
      for (scale in c(FALSE, TRUE)) {
        expected <- helpers$calcMarginalZScores(X, response, center = center, scale = scale)
        reset_counts()
        actual <- helpers$calcMarginalZScores(X, response, center = center, scale = scale,
                                             .x_stats = stats)
        expect_same(expected, actual, "Reuse must preserve every response shape and setting")
        stopifnot(counts$mean == 0L, counts$sd == 0L, counts$range == 0L)
      }
    }
  }

  # Reusing the wrong gene's values must fail these cases, even when matrix
  # dimensions match. Attribute changes also invalidate the input key.
  changed <- sweep(X, 2L, seq_len(ncol(X)) / 5 + 1, "*")
  changed <- sweep(changed, 2L, seq_len(ncol(X)) / 10, "+")
  reordered <- X[, rev(seq_len(ncol(X))), drop = FALSE]
  renamed <- X
  colnames(renamed) <- paste0("other", seq_len(ncol(X)))
  for (candidate in list(changed, reordered, renamed)) {
    tracked_X <- candidate
    expected <- helpers$computeMarginalZScores(candidate, y)
    reset_counts()
    actual <- helpers$computeMarginalZScores(candidate, y, .x_stats = stats)
    expect_same(expected, actual, "Stale or reordered input statistics must be rejected")
    stopifnot("Rejected input keys must use the original calculation" =
                counts$mean == 1L && counts$sd == 1L && counts$range == 1L)
  }
  tracked_X <- X
  invalid_stats <- list(
    list(), replace(stats, "n", list(nrow(X) + 1L)),
    replace(stats, "p", list(ncol(X) + 1L)),
    replace(stats, "finite", list(FALSE)),
    replace(stats, "x_mean", list(stats$x_mean[-1L])),
    replace(stats, "x_sd", list("invalid"))
  )
  for (invalid in invalid_stats) {
    reset_counts()
    expect_same(ordinary_z, helpers$computeMarginalZScores(X, y, .x_stats = invalid),
                "Malformed statistics must use the original calculation")
    stopifnot(counts$mean == 1L, counts$sd == 1L, counts$range == 1L)
  }

  # Statistics must not bypass native numerical fallbacks. Compare returned
  # values and conditions with the real installed implementation.
  capture <- function(operation, seed = 555L) {
    set.seed(seed)
    before <- .Random.seed
    notes <- character()
    output <- capture.output(value <- tryCatch(withCallingHandlers(operation(),
      warning = function(condition) {
        notes <<- c(notes, paste0("warning: ", conditionMessage(condition)))
        invokeRestart("muffleWarning")
      }, message = function(condition) {
        notes <<- c(notes, paste0("message: ", conditionMessage(condition)))
        invokeRestart("muffleMessage")
      }), error = function(condition) {
        structure(list(message = conditionMessage(condition), class = class(condition)),
                  class = "susier_test_error")
      }))
    list(value = value, notes = notes, output = output,
         rng_before = before, rng_after = .Random.seed)
  }
  constant <- X
  constant[, 2L] <- 0
  near_collinear <- X
  near_collinear[, 2L] <- X[, 1L] + seq_len(nrow(X)) * 1e-12
  fallback_cases <- list(
    constant = list(X = constant, y = y, scale = TRUE),
    offset = list(X = X + 1e8, y = y, scale = TRUE),
    tiny = list(X = X * 1e-170, y = y, scale = TRUE),
    perfect = list(X = X, y = X[, 1L], scale = TRUE),
    near_perfect = list(X = X, y = X[, 1L] + 1e-6 * y, scale = TRUE),
    constant_response = list(X = X, y = rep(1, nrow(X)), scale = TRUE),
    scale_ratio = list(X = X * 1e150, y = y * 1e-170, scale = FALSE)
  )
  native_z <- getFromNamespace("calc_z", "susieR")
  for (name in names(fallback_cases)) {
    case <- fallback_cases[[name]]
    expected <- capture(function() native_z(case$X, case$y, center = TRUE, scale = case$scale))
    reset_counts()
    actual <- capture(function() helpers$calcMarginalZScores(case$X, case$y,
      center = TRUE, scale = case$scale, .x_stats = make_stats(case$X)))
    expect_same(expected, actual, paste(name, "must preserve the native fallback"))
    stopifnot("The numerical guard must still invoke native calc_z" = counts$native_z == 1L)
  }

  # Record existing calculations in a real fit and check that the callback
  # receives only their exact raw results after the large cache is cleared.
  factory <- helpers$.createSusieFitEnvironment
  retained <- NULL
  fail_fit <- FALSE

  helpers$.createSusieFitEnvironment <- function(...) {
    private <- factory(...)
    retained <<- private
    if (is.null(private)) return(NULL)
    for (entry in list(c(".susier_native_colMeans", "fit_mean"),
                       c(".susier_native_compute_colSds", "fit_sd"))) {
      binding <- entry[[1L]]
      counter <- entry[[2L]]
      calculation <- private[[binding]]
      stopifnot("The fit must expose the audited native calculation" = is.function(calculation))
      private[[binding]] <- local({
        leaf <- calculation
        counter_name <- counter
        function(...) {
          counts[[counter_name]] <- counts[[counter_name]] + 1L
          leaf(...)
        }
      })
    }
    if (fail_fit) {
      saved <- private$.susier_native_get_ER2
      private$.susier_native_get_ER2 <- function(...) {
        saved(...)
        stop("Intentional fit failure after input statistics", call. = FALSE)
      }
    }
    private
  }
  reset_counts()
  captured_stats <- NULL
  direct_fit <- capture(function() helpers$runSusieWithFitOptimizations(X, y,
    L = 3L, max_iter = 30L, .postprocess = function(fit, transient) {
      stopifnot("The fit must finish before the callback starts" = counts$fit_sd > 0L,
                "The large cache must be empty before post-processing" =
                  length(ls(retained$.susier_fit_cache, all.names = TRUE)) == 0L)
      expect_same(stats, transient, "The callback must receive exact raw input statistics")
      captured_stats <<- transient
      fit$z <- helpers$computeMarginalZScores(X, y, .x_stats = transient)
      fit
    }))
  stopifnot(!inherits(direct_fit$value, "susier_test_error"),
            "A fit must compute its raw SD only once" = counts$fit_sd == 1L,
            "A fit must record its existing default mean calculation" = counts$fit_mean == 1L,
            "Post-fit scores must not repeat matrix scans" =
              counts$mean == 0L && counts$sd == 0L && counts$range == 1L,
            length(ls(retained$.susier_fit_cache, all.names = TRUE)) == 0L)
  expect_same(stats, captured_stats, "Copied statistics must survive private-cache cleanup")
  baseline_fit <- capture(function() {
    fit <- native(X, y, L = 3L, max_iter = 30L)
    fit$z <- helpers$computeMarginalZScores(X, y)
    fit
  })
  expect_same(baseline_fit, direct_fit, "Post-fit reuse must preserve the complete native fit")

  callback_error <- capture(function() helpers$runSusieWithFitOptimizations(X, y,
    L = 3L, max_iter = 30L, .postprocess = function(fit, transient) {
      stopifnot(length(ls(retained$.susier_fit_cache, all.names = TRUE)) == 0L)
      stop("Intentional post-fit failure", call. = FALSE)
    }))
  stopifnot(inherits(callback_error$value, "susier_test_error"),
            identical(callback_error$value$message, "Intentional post-fit failure"),
            length(ls(retained$.susier_fit_cache, all.names = TRUE)) == 0L)
  fail_fit <- TRUE
  callback_calls <- 0L
  fit_error <- capture(function() helpers$runSusieWithFitOptimizations(X, y,
    L = 3L, max_iter = 30L, .postprocess = function(fit, transient) {
      callback_calls <<- callback_calls + 1L
      fit
    }))
  stopifnot(inherits(fit_error$value, "susier_test_error"), callback_calls == 0L,
            identical(fit_error$value$message, "Intentional fit failure after input statistics"),
            length(ls(retained$.susier_fit_cache, all.names = TRUE)) == 0L)
  fail_fit <- FALSE

  # Native preprocessing can change the matrix. Its statistics must not be
  # used for the caller's original rows or columns.
  missing_y <- y
  missing_y[1L] <- NA_real_
  for (mode in list(list(y = y, extra = list(null_weight = .1)),
                    list(y = missing_y, extra = list(na.rm = TRUE)))) {
    callback_calls <- 0L
    arguments <- c(list(X = X, y = mode$y, L = 3L, max_iter = 30L), mode$extra)
    expected <- capture(function() do.call(native, arguments))
    arguments$.postprocess <- function(fit, transient) {
      callback_calls <<- callback_calls + 1L
      stopifnot("Changed native matrix dimensions must reject original-input reuse" =
                  is.null(transient),
                length(ls(retained$.susier_fit_cache, all.names = TRUE)) == 0L)
      fit
    }
    actual <- capture(function() do.call(helpers$runSusieWithFitOptimizations, arguments))
    expect_same(expected, actual, "Native row or null-column changes must preserve the fit")
    stopifnot(callback_calls == 1L)
  }

  # A source mismatch must complete the native fit before callback entry.
  fallback_helpers <- make_helpers()
  proxy <- new.env(parent = parent.env(namespace))
  list2env(namespace_values, envir = proxy)
  events <- character()
  proxy$susie <- function(...) {
    events <<- c(events, "fit-start")
    value <- native(...)
    events <<- c(events, "fit-end")
    value
  }
  fallback_helpers$asNamespace <- function(package) {
    if (identical(package, "susieR")) proxy else base::asNamespace(package)
  }
  fallback_fit <- capture(function() fallback_helpers$runSusieWithFitOptimizations(X, y,
    L = 3L, max_iter = 30L, .postprocess = function(fit, transient) {
      expect_same(c("fit-start", "fit-end"), events,
                  "Native fallback must finish before callback entry")
      stopifnot(is.null(transient))
      events <<- c(events, "callback")
      fit
    }))
  native_fit <- capture(function() native(X, y, L = 3L, max_iter = 30L))
  expect_same(native_fit, fallback_fit, "Source fallback must retain the complete native fit")
  expect_same(c("fit-start", "fit-end", "callback"), events,
              "Source fallback must run the callback exactly once")
  events <- character()
  expected_error <- capture(function() native(X, missing_y, L = 3L, max_iter = 30L))
  fallback_error <- capture(function() fallback_helpers$runSusieWithFitOptimizations(X,
    missing_y, L = 3L, max_iter = 30L, .postprocess = function(fit, transient) {
      events <<- c(events, "callback")
      fit
    }))
  expect_same(expected_error, fallback_error, "A failed native fallback must retain its error")
  expect_same("fit-start", events, "A failed native fit must not enter the callback")

  # Full pipeline fits must preserve serialized bytes, conditions and RNG.
  # The legacy path must retain native z-scores without post-fit replacement.
  fixtures <- list(
    ordinary = list(X = X, y = y, legacy = FALSE),
    constant = list(X = constant, y = cbind(y), legacy = FALSE),
    perfect = list(X = X, y = X[, 1L], legacy = FALSE),
    near_collinear = list(X = near_collinear, y = cbind(y), legacy = FALSE),
    legacy = list(X = near_collinear, y = cbind(y), legacy = TRUE)
  )
  for (name in names(fixtures)) {
    fixture <- fixtures[[name]]
    tracked_X <- fixture$X
    expected <- capture(function() {
      message("Fitting SuSiE model")
      fit <- native(fixture$X, fixture$y, L = 10,
        estimate_residual_variance = TRUE, estimate_prior_variance = TRUE,
        scaled_prior_variance = .1, verbose = TRUE,
        compute_univariate_zscore = fixture$legacy, min_abs_corr = .5)
      if (!fixture$legacy) {
        message("Calculating marginal z-scores")
        fit$z <- helpers$computeMarginalZScores(fixture$X, fixture$y)
      }
      fit$covariate_explained_variants <- helpers$emptyCovariateExplainedVariants()
      fit
    })
    reset_counts()
    actual <- capture(function() helpers$fitSusieModel(fixture$X, fixture$y,
                                                       legacy = fixture$legacy))
    expect_same(expected, actual, paste(name, "must preserve fit, z, conditions and RNG"))
    stopifnot("Pipeline fitting must compute raw input statistics once" =
                counts$fit_mean == 1L && counts$fit_sd == 1L,
              "The pipeline must consume transient statistics without another matrix scan" =
                counts$mean == 0L && counts$sd == 0L && counts$range == 1L,
              "The fixture must return a complete native fit" = inherits(actual$value, "susie"))
    expect_same(serialize(expected$value, NULL, version = 3),
                serialize(actual$value, NULL, version = 3),
                paste(name, "must preserve complete saved-model bytes"))
    stopifnot("The saved fit must not retain transient statistics" =
                !any(c(".x_stats", "x_stats", "X", ".susier_fit_cache") %in% names(actual$value)))
  }
  check_namespace()
  message("Marginal statistic reuse tests passed")
  invisible(TRUE)
}
