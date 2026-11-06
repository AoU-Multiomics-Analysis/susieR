# Exercise real native calculations. Spies count calls without changing values.
# Sourcing this file defines the suite but does not execute a fit.
run_susie_fit_cache_tests <- function(repo_dir = ".") {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
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
  serialized_same <- function(x, y) {
    identical(serialize(x, NULL, version = 3), serialize(y, NULL, version = 3))
  }
  functions_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH",
                              unset = file.path(repo_dir, "R", "utils"))
  helper_path <- file.path(functions_path, "SusiePerformance.R")
  helpers <- new.env(parent = globalenv())
  sys.source(helper_path, envir = helpers)
  expect_same(old_seed, get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
              "Sourcing fit helpers must preserve RNG state")
  required_api <- c("runSusieWithFitOptimizations", ".createSusieFitEnvironment",
                    ".susieFitSourceSupported")
  for (name in required_api) {
    if (!exists(name, envir = helpers, inherits = FALSE)) {
      stop("Missing fit optimization API: ", name, call. = FALSE)
    }
  }

  namespace <- asNamespace("susieR")
  native <- get("susie", envir = namespace, inherits = FALSE)
  namespace_names <- ls(namespace, all.names = TRUE)
  native_bindings <- mget(namespace_names, envir = namespace, inherits = FALSE)
  native_locks <- vapply(namespace_names, bindingIsLocked, logical(1), env = namespace)
  source_proxy <- function() {
    proxy <- new.env(parent = parent.env(namespace))
    list2env(native_bindings, envir = proxy)
    # Give the proxy its own version metadata. Do not change the real namespace.
    metadata <- new.env(parent = parent.env(proxy$.__NAMESPACE__.))
    list2env(as.list(proxy$.__NAMESPACE__., all.names = TRUE), envir = metadata)
    proxy$.__NAMESPACE__. <- metadata
    proxy
  }
  check_namespace <- function() {
    expect_same(namespace_names, ls(namespace, all.names = TRUE),
                "A fit must not add or remove native namespace bindings")
    expect_same(native_bindings,
                mget(namespace_names, envir = namespace, inherits = FALSE),
                "A fit must not change native namespace values")
    expect_same(native_locks,
                vapply(namespace_names, bindingIsLocked, logical(1), env = namespace),
                "A fit must preserve native binding locks")
    invisible(TRUE)
  }
  stopifnot("The audited local/image source must enable fit optimizations" =
              isTRUE(helpers$.susieFitSourceSupported(namespace)))
  private <- helpers$.createSusieFitEnvironment(namespace)
  stopifnot("The audited source must have a private fit environment" = is.environment(private),
            "Caches must be private to a fit" = is.environment(private$.susier_fit_cache))
  for (name in c("neg_loglik.individual", "gaussian_ser_lbf")) {
    expect_same(formals(get(name, namespace)), formals(get(name, private)),
                paste(name, "optimizer formals must remain unchanged"))
    expect_same(body(get(name, namespace)), body(get(name, private)),
                paste(name, "optimizer body must remain unchanged"))
  }

  # A real calculation must run on each miss. Returning the old result or
  # ignoring coefficient/uncertainty changes must fail these checks.
  calls <- new.env(parent = emptyenv())
  instrument <- function(environment, binding, counter) {
    calls[[counter]] <- 0L
    saved <- environment[[binding]]
    environment[[binding]] <- local({
      native_leaf <- saved
      counter_name <- counter
      function(...) {
        calls[[counter_name]] <- calls[[counter_name]] + 1L
        check_namespace()
        native_leaf(...)
      }
    })
  }
  native_xb <- get("compute_Xb", namespace)
  native_er2 <- get("get_ER2.individual", namespace)
  native_sds <- get("compute_colSds", namespace)
  native_stats <- get("compute_colstats", namespace)
  raw_x <- matrix(sin(seq_len(48L)), nrow = 12L,
                  dimnames = list(paste0("s", 1:12), paste0("v", 1:4)))
  stats <- native_stats(raw_x)
  x <- raw_x
  attr(x, "scaled:center") <- stats$cm
  attr(x, "scaled:scale") <- stats$csd
  attr(x, "d") <- stats$d
  b <- c(0.25, -0.5, 0.75, -0.125)
  instrument(private, ".susier_native_compute_Xb", "xb")
  expect_same(native_xb(x, b), private$compute_Xb(x, b), "Projection miss must use native output")
  expect_same(native_xb(x, b), private$compute_Xb(x, b), "Projection hit must use native output")
  stopifnot("Identical projection inputs must use one native calculation" = calls$xb == 1L)
  changed_b <- replace(b, 2L, b[2L] + 0.125)
  expect_same(native_xb(x, changed_b), private$compute_Xb(x, changed_b),
              "Changed coefficients must invalidate projection reuse")
  stopifnot(calls$xb == 2L)
  changed_x <- x
  changed_x[1L, 1L] <- changed_x[1L, 1L] + 0.5
  expect_same(native_xb(changed_x, b), private$compute_Xb(changed_x, b),
              "Changed matrix values must invalidate projection reuse")
  stopifnot(calls$xb == 3L)
  for (attribute in c("scaled:center", "scaled:scale")) {
    changed_x <- x
    attr(changed_x, attribute)[1L] <- attr(changed_x, attribute)[1L] + 0.25
    before <- calls$xb
    expect_same(native_xb(changed_x, b), private$compute_Xb(changed_x, b),
                paste(attribute, "must be part of the projection key"))
    stopifnot(calls$xb == before + 1L)
  }
  bounded <- helpers$.createSusieFitEnvironment(namespace, capacity = 2L)
  instrument(bounded, ".susier_native_compute_Xb", "bounded_xb")
  for (value in list(b, b + 0.1, b + 0.2, b)) {
    expect_same(native_xb(x, value), bounded$compute_Xb(x, value),
                "Bounded projection reuse must retain native values")
  }
  stopifnot("A projection cache must respect its capacity" = calls$bounded_xb == 4L)

  data <- structure(list(X = x, y = cos(seq_len(nrow(x)))), class = "individual")
  model <- list(alpha = rbind(c(.1, .2, .3, .4), c(.4, .3, .2, .1)),
                mu = rbind(b, -b), mu2 = rbind(b^2 + .1, b^2 + .2),
                Xr = native_xb(x, b / 2), X_theta = NULL, slot_weights = NULL,
                predictor_weights = stats$d)
  instrument(private, ".susier_native_get_ER2", "er2")
  expect_same(native_er2(data, model), private$get_ER2.individual(data, model),
              "ER2 miss must use the native matrix calculation")
  expect_same(native_er2(data, model), private$get_ER2.individual(data, model),
              "ER2 hit must preserve native output")
  stopifnot("Identical ER2 inputs must use one native calculation" = calls$er2 == 1L)
  changed_models <- list(
    mu2 = within(model, mu2 <- mu2 + 0.125),
    alpha = within(model, alpha[1L, 1L] <- alpha[1L, 1L] + 0.01),
    mu = within(model, mu[1L, 1L] <- mu[1L, 1L] + 0.01),
    Xr = within(model, Xr[1L] <- Xr[1L] + 0.01),
    X_theta = within(model, X_theta <- rep(0.1, nrow(data$X))),
    slot_weights = within(model, slot_weights <- c(0.9, 1.1)),
    predictor_weights = within(model, predictor_weights <- predictor_weights + 0.1)
  )
  stopifnot("The uncertainty fixture must change ER2 without changing posterior means" =
              !strict_same(native_er2(data, model), native_er2(data, changed_models$mu2)))
  for (name in names(changed_models)) {
    before <- calls$er2
    changed <- changed_models[[name]]
    expect_same(native_er2(data, changed), private$get_ER2.individual(data, changed),
                paste(name, "must be part of the ER2 key"))
    stopifnot(calls$er2 == before + 1L)
  }
  for (name in c("X", "y")) {
    changed <- data
    changed[[name]][1L] <- changed[[name]][1L] + 0.01
    before <- calls$er2
    expect_same(native_er2(changed, model), private$get_ER2.individual(changed, model),
                paste(name, "must be part of the ER2 key"))
    stopifnot(calls$er2 == before + 1L)
  }

  # A zero SD is raw data. Only the local standardization scale can replace it.
  constant_x <- cbind(raw_x, zero = 0, constant = 1)
  sd_private <- helpers$.createSusieFitEnvironment(namespace)
  instrument(sd_private, ".susier_native_compute_colSds", "sd")
  for (center in c(FALSE, TRUE)) {
    for (scale in c(FALSE, TRUE)) {
      expect_same(native_stats(constant_x, center, scale),
                  sd_private$compute_colstats(constant_x, center, scale),
                  paste("Column statistics must preserve native arithmetic", center, scale))
    }
  }
  expect_same(native_sds(constant_x), sd_private$compute_colSds(constant_x),
              "Cached raw SD values must retain zeros")
  stopifnot("Repeated column statistics must reuse one native SD scan" = calls$sd == 1L)

  # The constructor changes only a warning scan. It must preserve warning
  # indices, their order, and all constructor values at numerical boundaries.
  capture_constructor <- function(fun, matrix) {
    notes <- character()
    output <- capture.output(value <- tryCatch(withCallingHandlers(
      fun(matrix, cos(seq_len(nrow(matrix)))),
      warning = function(condition) {
        notes <<- c(notes, paste0("warning: ", conditionMessage(condition)))
        invokeRestart("muffleWarning")
      }, message = function(condition) {
        notes <<- c(notes, paste0("message: ", conditionMessage(condition)))
        invokeRestart("muffleMessage")
      }), error = function(condition) {
        list(error = conditionMessage(condition), class = class(condition))
      }))
    list(value = value, notes = notes, output = output)
  }
  constructor_cases <- list(
    constants = constant_x,
    large_offset = raw_x + 1e12,
    tiny = raw_x * 1e-200,
    huge = raw_x * 1e200,
    huge_constant = cbind(huge = rep(1e308, nrow(raw_x)), ordinary = raw_x[, 1L]),
    near_constant = cbind(1 + seq_len(12L) * .Machine$double.eps, raw_x),
    one_sample = raw_x[1L, , drop = FALSE],
    nonfinite = replace(raw_x, 1L, Inf)
  )
  for (name in names(constructor_cases)) {
    matrix <- constructor_cases[[name]]
    expect_same(capture_constructor(get("individual_data_constructor", namespace), matrix),
                capture_constructor(private$individual_data_constructor, matrix),
                paste(name, "constructor values and constant-column conditions must match native"))
  }

  # Body/API drift must use native fallback, even when the version is unchanged.
  guarded_names <- c("susie", "individual_data_constructor", "compute_colstats",
                     "compute_colSds", "compute_Xb", "compute_MXt",
                     "get_ER2.individual", "susie_workhorse")
  stopifnot("An unchanged source proxy must pass before mutation checks" =
              isTRUE(helpers$.susieFitSourceSupported(source_proxy())))
  metadata_proxy <- source_proxy()
  for (name in guarded_names) {
    metadata_proxy[[name]] <- eval(parse(
      text = deparse(get(name, namespace), control = "keepInteger"), keep.source = TRUE
    ), envir = namespace)
  }
  stopifnot("Source references must not disable an unchanged audited body" =
              isTRUE(helpers$.susieFitSourceSupported(metadata_proxy)))
  for (name in guarded_names) {
    proxy <- source_proxy()
    changed <- get(name, namespace)
    body(changed) <- substitute({ BODY; invisible(NULL) }, list(BODY = body(changed)))
    proxy[[name]] <- changed
    stopifnot("Changed audited source must disable private fits" =
                identical(helpers$.susieFitSourceSupported(proxy), FALSE),
              "An unsupported source must not create a private fit environment" =
                is.null(helpers$.createSusieFitEnvironment(proxy)))
  }
  proxy <- source_proxy()
  changed <- get("compute_Xb", namespace)
  formals(changed) <- c(formals(changed), alist(extra_argument = NULL))
  proxy$compute_Xb <- changed
  stopifnot("Changed audited APIs must disable private fits" =
              identical(helpers$.susieFitSourceSupported(proxy), FALSE))
  proxy <- source_proxy()
  proxy$.__NAMESPACE__.$spec["version"] <- "0.0.0"
  stopifnot("Unsupported package versions must disable private fits" =
              identical(helpers$.susieFitSourceSupported(proxy), FALSE),
            is.null(helpers$.createSusieFitEnvironment(proxy)))

  # Whole-fit fixtures preserve all fields, conditions, output, and RNG state.
  set.seed(20261106)
  n <- 80L
  p <- 40L
  ordinary <- matrix(rnorm(n * p), nrow = n,
                     dimnames = list(paste0("s", seq_len(n)), paste0("v", seq_len(p))))
  response <- .3 * ordinary[, 1L] - .2 * ordinary[, 20L] + rnorm(n)
  constants <- ordinary
  constants[, 4L] <- 0
  constants[, 6L] <- 1
  near <- ordinary
  near[, 2L] <- near[, 1L] + seq_len(n) * 1e-12
  near[, 3L] <- -near[, 1L]
  sampling <- ordinary
  sampling[, seq_len(10L)] <- ordinary[, rep(1L, 10L)]
  covariates <- cbind(binary = rep(c(0, 1), length.out = n), pc = rnorm(n), intercept = 1)
  raw_legacy <- ordinary
  raw_legacy[, 1L] <- covariates[, 1L]
  hat <- diag(n) - covariates %*% solve(crossprod(covariates)) %*% t(covariates)
  fixtures <- list(
    ordinary = list(X = ordinary, y = cbind(response)),
    constant = list(X = constants, y = cbind(response)),
    near_collinear = list(X = near, y = response),
    strong_signal = list(X = near,
                         y = cbind(3 * ordinary[, 1L] - 2 * ordinary[, 20L] + rnorm(n, sd = .3))),
    purity_subsampling = list(X = sampling,
                              y = cbind(3 * ordinary[, 1L] - 2 * ordinary[, 20L] + rnorm(n, sd = .3)),
                              n_purity = 2L),
    legacy_residuals = list(X = hat %*% raw_legacy, y = hat %*% response, legacy = TRUE)
  )
  capture_fit <- function(fun, fixture, extra = list()) {
    set.seed(12345)
    before <- .Random.seed
    notes <- character()
    arguments <- c(list(X = fixture$X, y = fixture$y, L = 10,
                        estimate_residual_variance = TRUE, estimate_prior_variance = TRUE,
                        scaled_prior_variance = .1, verbose = FALSE,
                        compute_univariate_zscore = isTRUE(fixture$legacy), min_abs_corr = .5,
                        n_purity = if (is.null(fixture$n_purity)) "auto" else fixture$n_purity), extra)
    output <- capture.output(value <- tryCatch(withCallingHandlers({
      fitted <- do.call(fun, arguments)
      if (!isTRUE(fixture$legacy)) {
        fitted$z <- helpers$computeMarginalZScores(fixture$X, fixture$y)
      }
      fitted
    }, warning = function(condition) {
      notes <<- c(notes, paste0("warning: ", conditionMessage(condition)))
      invokeRestart("muffleWarning")
    }, message = function(condition) {
      notes <<- c(notes, paste0("message: ", conditionMessage(condition)))
      invokeRestart("muffleMessage")
    }), error = function(condition) {
      structure(list(message = conditionMessage(condition), class = class(condition)),
                class = "susier_test_error")
    }))
    list(fit = value, rng_before = before, rng_after = .Random.seed,
         notes = notes, output = output)
  }
  compare_fit <- function(reference, candidate, label) {
    expect_same(reference, candidate, paste(label, "must preserve the complete fit and conditions"))
    stopifnot("Serialized fit bytes must match native" = serialized_same(reference$fit, candidate$fit))
    check_namespace()
  }
  references <- lapply(fixtures, function(fixture) capture_fit(native, fixture))
  stopifnot("The purity fixture must exercise native RNG use" =
              !strict_same(references$purity_subsampling$rng_before,
                           references$purity_subsampling$rng_after),
            "The strong-signal fixture must have credible sets" =
              length(references$strong_signal$fit$sets$cs) > 0L)
  for (name in names(fixtures)) {
    compare_fit(references[[name]],
                capture_fit(helpers$runSusieWithFitOptimizations, fixtures[[name]]), name)
  }

  # Count native leaf calls in a real fit. Parity alone must not allow a
  # permanently disabled cache or a clone that bypasses its private S3 methods.
  counted <- helpers$.createSusieFitEnvironment(namespace)
  requested <- new.env(parent = emptyenv())
  for (entry in list(c("compute_Xb", ".susier_native_compute_Xb", "fit_xb"),
                    c("get_ER2.individual", ".susier_native_get_ER2", "fit_er2"),
                    c("compute_colSds", ".susier_native_compute_colSds", "fit_sd"))) {
    leaf <- entry[1L]
    counter <- entry[3L]
    instrument(counted, entry[2L], counter)
    requested[[counter]] <- 0L
    saved <- counted[[leaf]]
    counted[[leaf]] <- local({
      calculation <- saved
      counter_name <- counter
      function(...) {
        requested[[counter_name]] <- requested[[counter_name]] + 1L
        calculation(...)
      }
    })
  }
  compare_fit(references$ordinary, capture_fit(counted$susie, fixtures$ordinary), "Private S3 fit")
  for (name in c("fit_xb", "fit_er2", "fit_sd")) {
    stopifnot("Each cache must reuse a real native result within a complete fit" =
                calls[[name]] > 0L && requested[[name]] > calls[[name]])
  }

  # Same-size fits with different genes/inputs must not share cache entries.
  second <- fixtures$ordinary
  second$X <- ordinary[, rev(seq_len(p)), drop = FALSE] * 1.25
  colnames(second$X) <- paste0("other", seq_len(p))
  second$y <- response + ordinary[, 8L]
  second_reference <- capture_fit(native, second)
  for (fixture in list(fixtures$ordinary, second, fixtures$ordinary)) {
    reference <- if (identical(fixture, second)) second_reference else references$ordinary
    compare_fit(reference, capture_fit(helpers$runSusieWithFitOptimizations, fixture),
                "Sequential A/B/A fit")
  }

  sparse <- fixtures$ordinary
  sparse$X <- Matrix::Matrix(sparse$X, sparse = TRUE)
  compare_fit(capture_fit(native, sparse), capture_fit(helpers$runSusieWithFitOptimizations, sparse),
              "Sparse native fallback")
  invalid <- fixtures$ordinary
  invalid$y[1L] <- NA_real_
  invalid_reference <- capture_fit(native, invalid)
  stopifnot(inherits(invalid_reference$fit, "susier_test_error"))
  compare_fit(invalid_reference, capture_fit(helpers$runSusieWithFitOptimizations, invalid),
              "Native input error")

  # Source drift must fall back to the real native entry point. The proxy
  # changes a leaf that fallback does not use; it never changes the namespace.
  fallback_helpers <- new.env(parent = globalenv())
  sys.source(helper_path, envir = fallback_helpers)
  proxy <- source_proxy()
  changed <- get("compute_Xb", namespace)
  body(changed) <- substitute({ BODY; invisible(NULL) }, list(BODY = body(changed)))
  proxy$compute_Xb <- changed
  fallback_helpers$asNamespace <- function(package) {
    if (identical(package, "susieR")) proxy else base::asNamespace(package)
  }
  compare_fit(references$ordinary,
              capture_fit(fallback_helpers$runSusieWithFitOptimizations, fixtures$ordinary),
              "Unsupported source fallback")

  # Capture a private environment only for this test so cleanup can be checked
  # after success and after a failure with populated projection/ER2 caches.
  scoped_helpers <- new.env(parent = globalenv())
  sys.source(helper_path, envir = scoped_helpers)
  factory <- scoped_helpers$.createSusieFitEnvironment
  retained <- NULL
  fail_after_cache <- FALSE
  scoped_helpers$.createSusieFitEnvironment <- function(...) {
    environment <- factory(...)
    retained <<- environment
    if (!is.null(environment) && fail_after_cache) {
      saved <- environment$.susier_native_get_ER2
      environment$.susier_native_get_ER2 <- function(...) {
        value <- saved(...)
        stop("Intentional failure after native ER2", call. = FALSE)
      }
    }
    environment
  }
  compare_fit(references$ordinary,
              capture_fit(scoped_helpers$runSusieWithFitOptimizations, fixtures$ordinary),
              "Private fit cleanup")
  stopifnot("Successful fits must clear retained input cache values" =
              length(ls(retained$.susier_fit_cache, all.names = TRUE)) == 0L)
  fail_after_cache <- TRUE
  failed <- capture_fit(scoped_helpers$runSusieWithFitOptimizations, fixtures$ordinary)
  stopifnot(inherits(failed$fit, "susier_test_error"),
            identical(failed$fit$message, "Intentional failure after native ER2"),
            "Failed fits must clear retained input cache values" =
              length(ls(retained$.susier_fit_cache, all.names = TRUE)) == 0L)
  check_namespace()
  fail_after_cache <- FALSE
  compare_fit(references$ordinary,
              capture_fit(scoped_helpers$runSusieWithFitOptimizations, fixtures$ordinary),
              "Valid fit after private failure")

  # Source must also leave an absent RNG seed absent.
  rm(".Random.seed", envir = .GlobalEnv)
  seedless_helpers <- new.env(parent = globalenv())
  sys.source(helper_path, envir = seedless_helpers)
  stopifnot("Sourcing fit helpers must not create an RNG seed" =
              !exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  check_namespace()
  message("SuSiE fit cache tests passed: native keys, exact fits, RNG, fallbacks, and cleanup")
  invisible(TRUE)
}
