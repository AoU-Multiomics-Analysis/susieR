# Benchmark prototype for susieR 0.16.6. Do not source this in production.
#
# calc_z() fits one regression with an intercept for each variant. For regular
# finite inputs, its t statistic can be calculated from centered sums instead.
# The residual degrees of freedom are n - 2, even when center = FALSE.
# Scaling a predictor by a positive value does not change this statistic.
#
# This prototype keeps the original function for inputs where QR rounding,
# rank detection, missing values, or singular regressions can affect behavior.
# It does not change SuSiE fitting, PIPs, credible sets, or purity calculations.

fast_calc_z <- function(X, Y, center = FALSE, scale = FALSE) {
  reference <- function(x, y) {
    getFromNamespace("calc_z", "susieR")(
      x, y, center = center, scale = scale
    )
  }

  if (!is.matrix(X) || !is.numeric(X) || ncol(X) < 1L ||
      length(center) != 1L || !is.logical(center) || is.na(center) ||
      length(scale) != 1L || !is.logical(scale) || is.na(scale) ||
      !requireNamespace("matrixStats", quietly = TRUE)) {
    return(reference(X, Y))
  }

  one_response <- function(y) {
    # Keep the original row-removal and dimension-dropping behavior for NA y.
    if (anyNA(y)) return(reference(X, y))
    x <- X
    n <- length(y)

    # range() checks the matrix without making an n x p logical array.
    if (n != nrow(x) || n <= 2L || !all(is.finite(range(x))) ||
        any(!is.finite(y))) {
      return(reference(x, y))
    }

    y_mean <- mean(y)
    yc <- y - y_mean
    syy <- sum(yc * yc)
    y_sd <- sqrt(syy / (n - 1))
    if (!is.finite(syy) || syy <= 0 ||
        abs(y_mean) >= 1e6 * y_sd) {
      return(reference(x, y))
    }

    x_mean <- colMeans(x)
    x_sd <- matrixStats::colSds(x)
    sxx <- (n - 1) * x_sd^2
    if (any(!is.finite(sxx)) || any(sxx <= 0) ||
        any(abs(x_mean) >= 1e6 * x_sd)) {
      return(reference(x, y))
    }

    # Do not allocate a centered n x p copy of x. Correct the product for
    # the small sum(yc) that can remain after floating-point centering.
    sxy <- as.numeric(crossprod(yc, x)) - sum(yc) * x_mean
    correlation <- (sxy / sqrt(sxx)) / sqrt(syy)
    residual_fraction <- 1 - correlation^2
    if (any(!is.finite(residual_fraction)) ||
        any(residual_fraction < 1e-8)) {
      return(reference(x, y))
    }

    unname(correlation * sqrt((n - 2) / residual_fraction))
  }

  if (is.null(dim(Y))) {
    return(one_response(Y))
  }
  if (!is.matrix(Y) || ncol(Y) < 1L) {
    return(reference(X, Y))
  }
  # calc_z() returns only a dim attribute; it does not retain input names.
  unname(do.call(cbind, lapply(seq_len(ncol(Y)), function(i) {
    one_response(Y[, i])
  })))
}

# The pipeline uses the default SuSiE intercept and standardization settings.
# A matrix response keeps the p x 1 z field used by extractResults().
fast_univariate_zscore <- function(X, y) {
  fast_calc_z(X, y, center = TRUE, scale = TRUE)
}

# Run small fixed fixtures against the installed implementation. No large
# matrix is made. Sourcing this file does not run the checks or change the RNG.
validate_fast_calc_z <- function(tolerance = 1e-10) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(20261106)
  x <- matrix(rnorm(48 * 7), nrow = 48,
              dimnames = list(paste0("s", seq_len(48)),
                              paste0("v", seq_len(7))))
  y <- 0.4 * x[, 1] + rnorm(48)
  fixtures <- list(
    finite = list(x = x, y = y),
    offsets = list(x = sweep(x, 2, seq_len(ncol(x)) * 3, "+"),
                   y = y + 8),
    mixed_scales = list(x = sweep(x, 2, 10^seq(-8, 8, length.out = 7), "*"),
                        y = y),
    response_matrix = list(x = x, y = cbind(first = y, second = y + x[, 2])),
    missing_y = list(x = x, y = replace(y, c(2, 7), NA_real_)),
    missing_y_one_column = list(x = x[, 1, drop = FALSE],
                                y = replace(y, c(2, 7), NA_real_)),
    constant_column = list(x = cbind(x, constant = 1), y = y),
    all_constant = list(x = matrix(1, 48, 3), y = y),
    constant_y = list(x = x, y = rep(1, 48)),
    zero_y = list(x = x, y = rep(0, 48)),
    large_offset = list(x = x + 1e10, y = y),
    perfect_correlation = list(x = x, y = x[, 1]),
    two_rows = list(x = x[1:2, , drop = FALSE], y = y[1:2]),
    one_row = list(x = x[1, , drop = FALSE], y = y[1]),
    missing_x = list(x = replace(x, 1, NA_real_), y = y)
  )

  capture <- function(fun) {
    notes <- character()
    value <- tryCatch(withCallingHandlers(fun(),
      warning = function(w) {
        notes <<- c(notes, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        notes <<- c(notes, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    ), error = function(e) e)
    list(value = value, notes = notes)
  }

  purrr::imap_dfr(fixtures, function(fixture, fixture_name) {
    purrr::map_dfr(c(FALSE, TRUE), function(center) {
      purrr::map_dfr(c(FALSE, TRUE), function(scale) {
        baseline <- capture(function() {
          getFromNamespace("calc_z", "susieR")(
            fixture$x, fixture$y, center = center, scale = scale
          )
        })
        candidate <- capture(function() {
          fast_calc_z(fixture$x, fixture$y, center = center, scale = scale)
        })
        error_baseline <- inherits(baseline$value, "error")
        error_candidate <- inherits(candidate$value, "error")
        both_numeric <- !error_baseline && !error_candidate
        same_attributes <- both_numeric &&
          identical(attributes(baseline$value), attributes(candidate$value))
        same_nonfinite <- both_numeric &&
          identical(is.na(baseline$value), is.na(candidate$value)) &&
          identical(is.nan(baseline$value), is.nan(candidate$value)) &&
          identical(is.infinite(baseline$value), is.infinite(candidate$value))
        finite <- if (both_numeric) is.finite(baseline$value) &
          is.finite(candidate$value) else logical()
        max_difference <- if (any(finite)) {
          max(abs(baseline$value[finite] - candidate$value[finite]))
        } else 0
        same_errors <- identical(error_baseline, error_candidate) &&
          (!error_baseline ||
             identical(conditionMessage(baseline$value),
                       conditionMessage(candidate$value)))
        tibble::tibble(
          fixture = fixture_name, center = center, scale = scale,
          max_abs_difference = max_difference,
          same_attributes = same_attributes,
          same_nonfinite = same_nonfinite,
          same_notes = identical(baseline$notes, candidate$notes),
          same_errors = same_errors,
          passed = same_errors && (error_baseline ||
            (same_attributes && same_nonfinite && max_difference <= tolerance))
        )
      })
    })
  })
}
