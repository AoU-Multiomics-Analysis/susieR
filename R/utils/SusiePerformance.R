# Calculate marginal regression t statistics without a centered copy of X.
# SuSiE fits each predictor with an intercept and n - 2 residual degrees of
# freedom. Positive predictor scaling does not change this statistic.
calcMarginalZScores <- function(X, Y, center = FALSE, scale = FALSE, .x_stats = NULL) {
  native <- getFromNamespace("calc_z", "susieR")
  reference <- function() native(X, Y, center = center, scale = scale)

  # Use the installed implementation for unsupported APIs and input types.
  # The image regression suite compares this helper with that implementation.
  if (!all(c("X", "Y", "center", "scale") %in% names(formals(native))) ||
      !is.matrix(X) || !is.numeric(X) || is.complex(X) || ncol(X) < 1L ||
      length(center) != 1L || !is.logical(center) || is.na(center) ||
      length(scale) != 1L || !is.logical(scale) || is.na(scale) ||
      !requireNamespace("matrixStats", quietly = TRUE)) {
    return(reference())
  }
  vector_response <- is.null(dim(Y))
  if (!is.numeric(Y) || is.complex(Y) ||
      (!vector_response && (!is.matrix(Y) || ncol(Y) < 1L))) {
    return(reference())
  }
  n <- nrow(X)
  response_rows <- if (vector_response) length(Y) else nrow(Y)
  # Statistics are transient output from this fit's private environment. Keep
  # the full matrix key so stale or reordered same-size inputs cannot reuse them.
  cached <- is.list(.x_stats) && identical(.x_stats$finite, TRUE) &&
    identical(.x_stats$n, n) && identical(.x_stats$p, ncol(X)) &&
    identical(.x_stats$X, X, num.eq = FALSE, single.NA = FALSE, attrib.as.set = FALSE) &&
    is.double(.x_stats$x_mean) && is.null(dim(.x_stats$x_mean)) &&
    length(.x_stats$x_mean) == ncol(X) &&
    is.double(.x_stats$x_sd) && is.null(dim(.x_stats$x_sd)) &&
    length(.x_stats$x_sd) == ncol(X)
  # range() avoids an n-by-p logical array when checking X.
  if (response_rows != n || n <= 2L || (!cached && !all(is.finite(range(X)))) ||
      any(!is.finite(Y))) {
    return(reference())
  }

  x_mean <- if (cached) .x_stats$x_mean else colMeans(X)
  x_sd <- if (cached) .x_stats$x_sd else matrixStats::colSds(X)
  sxx <- (n - 1) * x_sd^2
  minimum_sd <- sqrt(.Machine$double.xmin)
  # Keep native missing-data, rank-detection, and extreme-offset behavior.
  if (any(!is.finite(sxx)) || any(sxx <= 0) || any(x_sd < minimum_sd) ||
      any(abs(x_mean) >= 1e6 * x_sd)) {
    return(reference())
  }

  responses <- if (vector_response) matrix(Y, ncol = 1L) else Y
  y_mean <- colMeans(responses)
  yc <- sweep(responses, 2, y_mean, "-")
  syy <- colSums(yc * yc)
  y_sd <- sqrt(syy / (n - 1))
  if (any(!is.finite(syy)) || any(syy <= 0) || any(y_sd < minimum_sd) ||
      any(abs(y_mean) >= 1e6 * y_sd)) {
    return(reference())
  }
  # Native standard errors can underflow or overflow at extreme scale ratios.
  # Keep that behavior, even when the scale-free statistic would be finite.
  if (!scale &&
      (min(y_sd) / max(x_sd) < minimum_sd * sqrt(n) ||
       max(y_sd) / min(x_sd) > sqrt(.Machine$double.xmax) / sqrt(n))) {
    return(reference())
  }

  # Correct the small sum(yc) that can remain after floating-point centering.
  sxy <- crossprod(X, yc) - outer(x_mean, colSums(yc))
  correlation <- sweep(sweep(sxy, 1, sqrt(sxx), "/"), 2, sqrt(syy), "/")
  residual_fraction <- 1 - correlation^2
  if (any(!is.finite(residual_fraction)) || any(residual_fraction < 1e-8)) {
    return(reference())
  }
  z <- unname(correlation * sqrt((n - 2) / residual_fraction))
  if (vector_response) as.vector(z) else z
}

# Match the pipeline's intercept and predictor-standardization settings.
computeMarginalZScores <- function(X, y, .x_stats = NULL) {
  calcMarginalZScores(X, y, center = TRUE, scale = TRUE, .x_stats = .x_stats)
}

# These source signatures identify the audited individual-data implementation
# in the pinned SuSiE image. A version number alone does not identify its code.
.susieFitSourceSupported <- function(ns) {
  tryCatch({
    if (!is.environment(ns) ||
        !identical(as.character(getNamespaceVersion(ns)), "0.16.6") ||
        !requireNamespace("digest", quietly = TRUE) ||
        !requireNamespace("matrixStats", quietly = TRUE)) return(FALSE)
    signatures <- c(
      susie = "e784df68f623b4c0f3785fb4faaeea5bfa179531d21f9a1f578acfb16ac142e9",
      individual_data_constructor = "a5f82fdeb850fd7adfb94e00a8a0f64aa96b78d04c55003c68925ed45631e893",
      compute_colstats = "496859f9c1759081e623d6642c810f5e648e7cc1aa94382cec928e0300eb9537",
      compute_colSds = "8a0061c5b201e5934b78ec10a64d2c1f6226563718ae9592eeee5fc898b87f81",
      compute_Xb = "e80b9ec503ce1aa42971878da7d117199bbb51779a172086bd6c2cc3853f456f",
      get_ER2.individual = "8914694f862af8b6b6ba2c561c3b2b9bee6760a7c75b301188cb1c1283bf1b7a",
      susie_workhorse = "69b716e0f96afd9341bdaaf02608553430d312db39739b2e02ceb1b4f80038ac",
      compute_MXt = "64c9f659fe153904ad524dc30b4b02f154b5748222b2b377c7be75ce3c72d128"
    )
    all(vapply(names(signatures), function(name) {
      fun <- get0(name, envir = ns, inherits = FALSE)
      if (!is.function(fun) || is.primitive(fun)) return(FALSE)
      # Parsing without source references also strips metadata from nested
      # anonymous functions. Compare code, not bytecode or source locations.
      expression <- parse(text = deparse(utils::removeSource(fun), width.cutoff = 500L),
                          keep.source = FALSE)[[1L]]
      if (!is.call(expression) || !identical(expression[[1L]], quote(`function`))) {
        return(FALSE)
      }
      actual <- digest::digest(list(expression[[2L]], expression[[3L]]),
                               algo = "sha256", serializeVersion = 2)
      identical(actual, unname(signatures[[name]]))
    }, logical(1)))
  }, error = function(condition) FALSE)
}

# Rebind the complete private closure graph so its S3 callers use private
# methods too. No binding or installed file in the native namespace is changed.
.createSusieFitEnvironment <- function(ns = asNamespace("susieR"), capacity = 22L) {
  if (!.susieFitSourceSupported(ns) || length(capacity) != 1L ||
      !is.numeric(capacity) || !is.finite(capacity) || capacity < 1 ||
      capacity > .Machine$integer.max || capacity != floor(capacity)) return(NULL)
  private <- new.env(parent = ns)
  for (name in ls(ns, all.names = TRUE)) {
    fun <- get(name, envir = ns, inherits = FALSE)
    if (is.function(fun) && identical(environment(fun), ns)) {
      environment(fun) <- private
      assign(name, fun, envir = private)
    }
  }
  cache <- new.env(parent = emptyenv())
  cache$capacity <- as.integer(capacity)
  cache$projections <- list()
  private$.susier_fit_cache <- cache
  private$.susier_same <- function(x, y) {
    identical(x, y, num.eq = FALSE, single.NA = FALSE, attrib.as.set = FALSE)
  }
  private$.susier_native_compute_Xb <- private$compute_Xb
  private$.susier_native_get_ER2 <- private$get_ER2.individual
  private$.susier_native_compute_colSds <- private$compute_colSds
  private$.susier_native_colMeans <- get("colMeans", envir = ns, inherits = TRUE)

  projection <- function(X, b) {
    cache <- .susier_fit_cache
    if (is.null(cache$projection_X) || !.susier_same(X, cache$projection_X)) {
      cache$projection_X <- X
      cache$projections <- list()
    }
    for (i in seq_along(cache$projections)) {
      entry <- cache$projections[[i]]
      if (.susier_same(b, entry$b)) {
        cache$projections <- c(cache$projections[-i], list(entry))
        return(entry$value)
      }
    }
    value <- .susier_native_compute_Xb(X, b)
    cache$projections <- tail(c(cache$projections, list(list(b = b, value = value))),
                              cache$capacity)
    value
  }
  environment(projection) <- private
  private$compute_Xb <- projection

  residual_error <- function(data, model) {
    cache <- .susier_fit_cache
    key <- list(data$X, data$y, model$alpha, model$mu, model$mu2, model$Xr,
                model$X_theta, model$slot_weights, model$predictor_weights)
    if (!is.null(cache$er2_key) && .susier_same(key, cache$er2_key)) {
      return(cache$er2_value)
    }
    value <- .susier_native_get_ER2(data, model)
    cache$er2_key <- key
    cache$er2_value <- value
    value
  }
  environment(residual_error) <- private
  private$get_ER2.individual <- residual_error

  column_sds <- function(X) {
    cache <- .susier_fit_cache
    if (!is.null(cache$sd_X) && .susier_same(X, cache$sd_X)) return(cache$sd_value)
    value <- .susier_native_compute_colSds(X)
    cache$sd_X <- X
    # Keep the raw zeros. Native compute_colstats changes only its local copy
    # when it replaces zero SDs with one for predictor standardization.
    cache$sd_value <- value
    value
  }
  environment(column_sds) <- private
  private$compute_colSds <- column_sds

  warning_variances <- function(X) {
    values <- matrixStats::colVars(X)
    exceptional <- which(!is.finite(values))
    if (length(exceptional)) {
      # Overflow can give NaN here but Inf in native var(), which changes the
      # warning rule. Recheck only these columns without an apply-sized copy.
      values[exceptional] <- vapply(exceptional, function(i) var(X[, i]), numeric(1))
    }
    values
  }
  environment(warning_variances) <- private
  private$.susier_warning_variances <- warning_variances
  replacements <- 0L
  replace_warning_scan <- function(expression) {
    if (identical(expression, quote(apply(X, 2, var)))) {
      replacements <<- replacements + 1L
      return(quote(.susier_warning_variances(X)))
    }
    if (is.call(expression)) {
      for (i in seq_along(expression)) {
        expression[[i]] <- replace_warning_scan(expression[[i]])
      }
    }
    expression
  }
  constructor <- private$individual_data_constructor
  body(constructor) <- replace_warning_scan(body(constructor))
  if (replacements != 1L) return(NULL)
  private$individual_data_constructor <- constructor

  record_means <- function(X) {
    value <- .susier_native_colMeans(X)
    .susier_fit_cache$mean_X <- X
    .susier_fit_cache$mean_value <- value
    value
  }
  environment(record_means) <- private
  private$.susier_record_means <- record_means
  replacements <- 0L
  replace_default_mean <- function(expression) {
    if (identical(expression, quote(colMeans(X)))) {
      replacements <<- replacements + 1L
      return(quote(.susier_record_means(X)))
    }
    if (is.call(expression)) {
      for (i in seq_along(expression)) {
        expression[[i]] <- replace_default_mean(expression[[i]])
      }
    }
    expression
  }
  column_stats <- private$compute_colstats
  body(column_stats) <- replace_default_mean(body(column_stats))
  if (replacements != 1L) return(NULL)
  private$compute_colstats <- column_stats
  private
}

# Use caches only for the audited dense, ordinary single-effect model. All
# native arguments and optimizer/convergence settings pass through unchanged.
runSusieWithFitOptimizations <- function(X, y, ..., .postprocess = NULL) {
  ns <- asNamespace("susieR")
  native <- get("susie", envir = ns, inherits = FALSE)
  if (!is.null(.postprocess) && !is.function(.postprocess)) {
    stop(".postprocess must be NULL or a function", call. = FALSE)
  }
  finish <- function(fitted, stats = NULL) {
    if (is.null(.postprocess)) fitted else .postprocess(fitted, stats)
  }
  use_native <- function() {
    # Complete the fit before invoking a callback, including native fallbacks.
    fitted <- native(X, y, ...)
    finish(fitted)
  }
  supported <- is.matrix(X) && is.double(X) && !is.object(X) &&
    nrow(X) > 1L && ncol(X) > 0L && is.null(attr(X, "matrix.type")) &&
    all(is.finite(range(X)))
  if (!supported) return(use_native())
  arguments <- list(...)
  defaults <- list(unmappable_effects = "none", estimate_prior_method = "optim",
                   estimate_residual_method = "MoM", refine = FALSE,
                   slot_prior = NULL, L_greedy = NULL, model_init = NULL, s_init = NULL)
  for (name in names(defaults)) {
    if (name %in% names(arguments) && !identical(arguments[[name]], defaults[[name]])) {
      return(use_native())
    }
  }
  private <- .createSusieFitEnvironment(ns)
  if (is.null(private)) return(use_native())
  cache <- private$.susier_fit_cache
  on.exit(rm(list = ls(cache, all.names = TRUE), envir = cache), add = TRUE)
  fitted <- private$susie(X, y, ...)
  stats <- NULL
  if (!is.null(.postprocess) &&
      private$.susier_same(cache$mean_X, X) && private$.susier_same(cache$sd_X, X)) {
    stats <- list(X = X, n = nrow(X), p = ncol(X), finite = TRUE,
                  x_mean = cache$mean_value, x_sd = cache$sd_value)
  }
  # The callback needs only raw input statistics. Release scaled matrices and
  # effect caches now; on.exit also covers fitting and callback failures.
  rm(list = ls(cache, all.names = TRUE), envir = cache)
  finish(fitted, stats)
}
