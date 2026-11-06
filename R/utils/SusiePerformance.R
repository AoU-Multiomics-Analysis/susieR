# Calculate marginal regression t statistics without a centered copy of X.
# SuSiE fits each predictor with an intercept and n - 2 residual degrees of
# freedom. Positive predictor scaling does not change this statistic.
calcMarginalZScores <- function(X, Y, center = FALSE, scale = FALSE) {
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
  # range() avoids an n-by-p logical array when checking X.
  if (response_rows != n || n <= 2L || !all(is.finite(range(X))) ||
      any(!is.finite(Y))) {
    return(reference())
  }

  x_mean <- colMeans(X)
  x_sd <- matrixStats::colSds(X)
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
computeMarginalZScores <- function(X, y) {
  calcMarginalZScores(X, y, center = TRUE, scale = TRUE)
}
