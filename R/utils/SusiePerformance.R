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

# The SS path uses only the audited ordinary model and exported native helpers.
# A different installed implementation keeps the individual-data fit path.
.susieSufficientStatsSourceSupported <- function(ns = asNamespace("susieR")) {
  if (!.susieFitSourceSupported(ns)) return(FALSE)
  signatures <- c(
    susie_ss = "1dc18a6450b4712b7c1a0657f0aee10c2007c9205444c88a34e1f35ac9dadc97",
    sufficient_stats_constructor = "6b6a8c363496c853d6902f6c105d3725a19ae7e64bafc208f1fe207aa60c8e79",
    susie_get_cs = "514ef2bb8c067d227b80d0959744a01529e15e104340acf53b954b6552113508",
    get_fitted.individual = "10d86e9c0af42a4285662c95a397a4ffe0057d9350de13781179f4ed8faa9360",
    calculate_posterior_moments.ss = "d94182619e1cfe62972f871b3a643f01f12b659b5ea6fe9491a18bf4b74bdb6f",
    cleanup_model.ss = "ae6fa6fa8fdc481f7d0bd7ede48a72cf2306768ca0adc81e9df228ee3abb0ac8",
    compute_kl.ss = "fce7bc195f39a843f63e6f475e90f974f2f1b19c0f0892bf33a43ee38260ce4f",
    compute_residuals.ss = "9a713343d8bde116b558611868e565169723919f512a55829339cfb0354422ad",
    compute_ser_statistics.ss = "8c7e2340a60e96fc35e95ca873c3ef42d0b3244a051dce830361a930279e9e0a",
    configure_data.ss = "ab7864f6dc28564cc857000f436c5f003d353dde2f833990d9d7d448ebdd3333",
    Eloglik.ss = "d027639845016dc073759a06cf63337e4726de28e557ac38fc8a72dd5f0dfa02",
    get_cs.ss = "9a0278f1ec64285e7f58310eef23251eaa0a44f899faf084f2320f04a47e17b0",
    get_ER2.ss = "77893ac8d2e104052c5b6216dfb96d1cebc71b7333e14b370fb9a26c549f56a4",
    get_fitted.ss = "f893ecfffa4c50e899ddd78538bafcb4d7a48627a0e7c467f3ab4871601ae259",
    get_intercept.ss = "18a46dd65a512a94b7fc4f4c6a4c1f6d51ecc8444a77c83cc8acfee19fce1f6c",
    get_scale_factors.ss = "8f82efb8437504bcae1ab87ed87c8817c1736f0899f4a970f9bb96ba0ff4cf8e",
    get_var_y.ss = "0e85b3a0589ad17d6a017130e6210465665740e2081f085fcb4ec48b5649fdfa",
    get_variable_names.ss = "6e4958d22357d6701c128fcbc9f479d2bf56fa89bc1263083a432405b25e1491",
    get_zscore.ss = "55735c9c08f20da36bbf5e9f1159abfce473f7fc5ab3711802e247540cd5f2a3",
    initialize_fitted.ss = "30081faddad32dff98a3b6cc3e8af4546196de95e4e01dab70e3347e4769ecf1",
    initialize_susie_model.ss = "96b2700f44e30437a80c2dd8ea2a8b4128d2000bcab2a23d33c1c21cc878e177",
    loglik.ss = "bae822a65293a1e09cd097e155c1f48447df84453726af214542bdd0e14055d6",
    neg_loglik.ss = "d951344cba18a2c28efa950a73a55e611baeb7111885706915eacea80b7c18b4",
    SER_posterior_e_loglik.ss = "9f3b4de94791331c45e4d5658c253a7794848574a0f6a1def739447e4430b224",
    track_ibss_fit.ss = "ec2ae846914b9b80a8d58c99a7439754d3edcf149fd32bc02bb1abbe85f4dbf8",
    update_derived_quantities.ss = "1f41a000f5aacfdc1747bca8b32004a1d3d7fbb6cd017aedfb1c4e8a44c1efe6",
    update_fitted_values.ss = "de8813e2e4d92459101c62dd3a7f578843f336afd3bb5f107fbd45906dc6a6db",
    update_variance_components.ss = "3ef2a8dd3e709dc072e24e6a896adaedbc40b6891aa5f92be37432d9179e90b5",
    validate_prior.ss = "b0b097240ccc1af2781ab3071da5f0ccedc495a7a10b366c9f1493c46101cf76",
    assign_names = "269d5ab36c3c53dbd83fab261ff9cc2c8a570b607ff2474787e9e03e84225607")
  tryCatch(all(vapply(names(signatures), function(name) {
    fun <- get0(name, envir = ns, inherits = FALSE)
    if (!is.function(fun) || is.primitive(fun)) return(FALSE)
    expression <- parse(text = deparse(utils::removeSource(fun), width.cutoff = 500L),
                        keep.source = FALSE)[[1L]]
    actual <- digest::digest(list(expression[[2L]], expression[[3L]]),
                             algo = "sha256", serializeVersion = 2)
    identical(actual, unname(signatures[[name]]))
  }, logical(1))), error = function(condition) FALSE)
}

finalizeSusieSufficientStatistics <- function(fitted, X, y, statistics) {
  # Native SS purity uses full Xcorr. Use the individual-data API instead so
  # large sets use the same default purity sampling and consume the same RNG.
  fitted$sets <- susieR::susie_get_cs(fitted, X = X, coverage = 0.95,
    min_abs_corr = 0.5, n_purity = "auto")
  # Preserve saved-fit predictions. Use all effects, as native fitted values do;
  # posterior_mean() deliberately removes effects below its prior tolerance.
  beta <- colSums(fitted$alpha * fitted$mu) / fitted$X_column_scale_factors
  fitted$Xr <- as.vector(X %*% beta - sum(statistics$X_colmeans * beta))
  fitted$fitted <- fitted$Xr + mean(y)
  names(fitted$fitted) <- if (is.null(names(drop(y)))) rownames(X) else names(drop(y))
  fitted$XtXr <- NULL
  fitted
}

runSusieWithSufficientStatistics <- function(X, y, sufficient_statistics) {
  if (!.susieSufficientStatsSourceSupported(asNamespace("susieR"))) {
    message("USE individual-data SuSiE: sufficient-statistics source is not supported")
    return(NULL)
  }
  statistics <- sufficient_statistics
  if (!is.list(statistics) || !is.numeric(statistics$n) ||
      length(statistics$n) != 1L || is.na(statistics$n) ||
      statistics$n != nrow(X) || length(y) != nrow(X)) {
    stop("Sufficient statistics sample count does not match the fit", call. = FALSE)
  }
  if (!is.matrix(statistics$XtX) || !is.double(statistics$XtX) ||
      !identical(dim(statistics$XtX), c(ncol(X), ncol(X))) ||
      !is.numeric(statistics$X_colmeans) || length(statistics$X_colmeans) != ncol(X)) {
    stop("Sufficient statistics matrix dimension does not match the fit", call. = FALSE)
  }
  if (!identical(colnames(statistics$XtX), colnames(X)) ||
      !identical(rownames(statistics$XtX), colnames(X))) {
    stop("Sufficient statistics variant order does not match the fit", call. = FALSE)
  }
  if (!all(is.finite(statistics$XtX)) || any(!is.finite(statistics$X_colmeans)) ||
      any(!is.finite(y)) || any(diag(statistics$XtX) <= 0)) {
    message("USE individual-data SuSiE: sufficient statistics contain an unusable column")
    return(NULL)
  }
  yc <- drop(y) - mean(y)
  message("Fitting SuSiE model with reused XtX and susie_ss")
  fitted <- susieR::susie_ss(XtX = statistics$XtX,
    Xty = drop(yc %*% X), yty = sum(yc^2), n = nrow(X),
    X_colmeans = statistics$X_colmeans, y_mean = mean(y),
    L = 10, standardize = TRUE, scaled_prior_variance = 0.1,
    estimate_residual_variance = TRUE, estimate_prior_variance = TRUE,
    residual_variance_lowerbound = stats::var(drop(y)) / 10000,
    coverage = NULL, min_abs_corr = 0.5, verbose = TRUE)
  finalizeSusieSufficientStatistics(fitted, X, y, statistics)
}
