# Compare the bulk statistic with the installed native implementation.
# Wrong centering, residual degrees of freedom, shape, or fallback handling
# must fail these tests. Sourcing this file does not execute the tests.
run_susie_performance_tests <- function(repo_dir = ".") {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  helpers <- new.env(parent = globalenv())
  functions_dir <- Sys.getenv("SUSIER_FUNCTIONS_PATH",
                             unset = file.path(repo_dir, "R", "utils"))
  helper_path <- file.path(functions_dir, "SusiePerformance.R")
  stopifnot("Marginal-z helper must exist" = file.exists(helper_path))
  source_seed <- if (had_seed) old_seed else NULL
  sys.source(helper_path, envir = helpers)
  stopifnot("Sourcing must preserve the RNG" = identical(
    source_seed, get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  ))

  native <- getFromNamespace("calc_z", "susieR")
  capture <- function(fun) {
    notes <- character()
    value <- tryCatch(withCallingHandlers(fun(),
      warning = function(w) {
        notes <<- c(notes, paste0("warning: ", conditionMessage(w)))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        notes <<- c(notes, paste0("message: ", conditionMessage(m)))
        invokeRestart("muffleMessage")
      }
    ), error = function(e) e)
    list(value = value, notes = notes)
  }
  compare <- function(baseline, candidate, label, tolerance = 1e-10) {
    fail <- function(detail) stop(label, ": ", detail, call. = FALSE)
    if (!identical(baseline$notes, candidate$notes)) fail("different warnings or messages")
    baseline_error <- inherits(baseline$value, "error")
    candidate_error <- inherits(candidate$value, "error")
    if (!identical(baseline_error, candidate_error)) fail("different error behavior")
    if (baseline_error) {
      if (!identical(conditionMessage(baseline$value),
                     conditionMessage(candidate$value))) fail("different error message")
      return(invisible(0))
    }
    a <- baseline$value
    b <- candidate$value
    if (!identical(typeof(a), typeof(b))) fail("different output type")
    if (!identical(attributes(a), attributes(b))) fail("different output attributes")
    if (!identical(length(a), length(b))) fail("different output length")
    if (!identical(is.na(a), is.na(b)) || !identical(is.nan(a), is.nan(b)) ||
        !identical(is.infinite(a), is.infinite(b))) fail("different nonfinite masks")
    finite <- is.finite(a) & is.finite(b)
    difference <- if (any(finite)) max(abs(a[finite] - b[finite])) else 0
    if (difference > tolerance) fail(paste("maximum difference", difference))
    if (!identical(a[is.infinite(a)], b[is.infinite(b)])) fail("different infinite signs")
    invisible(difference)
  }

  set.seed(20261106)
  x <- matrix(rnorm(48 * 7), nrow = 48,
              dimnames = list(paste0("s", seq_len(48)), paste0("v", seq_len(7))))
  y <- 0.4 * x[, 1] + rnorm(48)
  names(y) <- rownames(x)
  fixtures <- list(
    finite = list(x = x, y = y),
    offsets = list(x = sweep(x, 2, seq_len(ncol(x)) * 3, "+"), y = y + 8),
    mixed_scales = list(x = sweep(x, 2, 10^seq(-8, 8, length.out = 7), "*"), y = y),
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
    missing_x = list(x = replace(x, 1, NA_real_), y = y),
    one_column = list(x = x[, 1, drop = FALSE], y = y),
    one_response_column = list(x = x, y = cbind(trait = y)),
    integer_matrix = list(x = matrix(seq_len(48 * 7) %% 5L, 48), y = y),
    near_perfect = list(x = x, y = x[, 1] + 1e-6 * y),
    huge_x = list(x = x * 1e200, y = y),
    tiny_x = list(x = x * 1e-200, y = y),
    subnormal_x = list(x = x * 1e-160, y = y),
    near_subnormal_x = list(x = x * 1e-155, y = y),
    huge_y = list(x = x, y = y * 1e200),
    tiny_y = list(x = x, y = y * 1e-200),
    subnormal_y = list(x = x, y = y * 1e-160),
    large_scale_ratio = list(x = x * 1e-100, y = y * 1e100),
    small_scale_ratio = list(x = x * 1e100, y = y * 1e-100),
    large_y_offset = list(x = x, y = y + 1e12),
    infinite_x = list(x = replace(x, 1, Inf), y = y),
    infinite_y = list(x = x, y = replace(y, 1, Inf)),
    wrong_rows = list(x = x, y = y[-1]),
    multiple_y_with_missing = list(x = x, y = cbind(y, replace(y, 1, NA_real_))),
    sparse = list(x = Matrix::Matrix(x, sparse = TRUE), y = y)
  )
  differences <- numeric()
  for (fixture_name in names(fixtures)) {
    fixture <- fixtures[[fixture_name]]
    for (center in c(FALSE, TRUE)) {
      for (scale in c(FALSE, TRUE)) {
        baseline <- capture(function() native(fixture$x, fixture$y,
                                               center = center, scale = scale))
        candidate <- capture(function() helpers$calcMarginalZScores(
          fixture$x, fixture$y, center = center, scale = scale
        ))
        differences <- c(differences, compare(baseline, candidate,
          paste(fixture_name, center, scale)))
      }
    }
    compare(capture(function() native(fixture$x, fixture$y, center = TRUE, scale = TRUE)),
            capture(function() helpers$computeMarginalZScores(fixture$x, fixture$y)),
            paste(fixture_name, "pipeline wrapper"))
  }

  # Native calc_z drops input names; keep its p-vector and p-by-q matrix shapes.
  z_vector <- helpers$computeMarginalZScores(x, y)
  z_matrix <- helpers$computeMarginalZScores(x, cbind(trait = y))
  stopifnot(is.double(z_vector), is.null(attributes(z_vector)),
            length(z_vector) == ncol(x),
            identical(dim(z_matrix), c(ncol(x), 1L)),
            is.null(dimnames(z_matrix)), all(is.finite(z_matrix)))

  # The helper must not create or advance the RNG, including fallback paths.
  rng_seed <- get(".Random.seed", envir = .GlobalEnv)
  invisible(helpers$computeMarginalZScores(x, y))
  invisible(capture(function() helpers$computeMarginalZScores(x, rep(1, nrow(x)))))
  stopifnot(identical(rng_seed, get(".Random.seed", envir = .GlobalEnv)))
  rm(".Random.seed", envir = .GlobalEnv)
  invisible(helpers$computeMarginalZScores(x, y))
  stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  assign(".Random.seed", rng_seed, envir = .GlobalEnv)

  # Compute the true native fit as well as its directly calculated marginal z.
  fit <- susieR::susie(x, cbind(y), L = 10, compute_univariate_zscore = TRUE)
  compare(list(value = fit$z, notes = character()),
          capture(function() helpers$computeMarginalZScores(x, cbind(y))),
          "complete native fit")
  message("Marginal-z tests passed: ", length(differences),
          " generic fixtures; maximum absolute difference = ",
          format(max(differences), scientific = TRUE))
  invisible(TRUE)
}
