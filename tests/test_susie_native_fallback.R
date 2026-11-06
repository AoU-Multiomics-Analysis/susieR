# A changed source guard must use real native fits without caches or input stats.
run_susie_native_fallback_tests <- function(repo_dir = ".") {
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
  helpers <- new.env(parent = globalenv())
  functions_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH",
                               unset = file.path(repo_dir, "R", "utils"))
  sys.source(file.path(functions_path, "SusiePerformance.R"), envir = helpers)
  namespace <- asNamespace("susieR")
  native <- get("susie", envir = namespace, inherits = FALSE)
  native_names <- ls(namespace, all.names = TRUE)
  native_values <- mget(native_names, envir = namespace, inherits = FALSE)
  native_locks <- vapply(native_names, bindingIsLocked, logical(1), env = namespace)
  check_namespace <- function() {
    stopifnot(strict_same(native_names, ls(namespace, all.names = TRUE)),
              strict_same(native_values, mget(native_names, envir = namespace, inherits = FALSE)),
              strict_same(native_locks,
                          vapply(native_names, bindingIsLocked, logical(1), env = namespace)))
  }
  # Override only the helper's guard. The installed namespace stays untouched.
  helpers$.susieFitSourceSupported <- function(ns) FALSE
  factory <- helpers$.createSusieFitEnvironment
  helpers$.createSusieFitEnvironment <- function(...) {
    private <- factory(...)
    if (!is.null(private)) stop("Unsupported source must not create fit caches", call. = FALSE)
    private
  }
  capture_fit <- function(fun, fixture, callback = NULL) {
    set.seed(12345)
    before <- .Random.seed
    notes <- character()
    arguments <- list(X = fixture$X, y = fixture$y, L = 3L, max_iter = 30L,
                      verbose = FALSE, min_abs_corr = .5, n_purity = 2L)
    if (!is.null(callback)) arguments$.postprocess <- callback
    output <- capture.output(value <- tryCatch(withCallingHandlers(
      do.call(fun, arguments),
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
  set.seed(20261106)
  X <- matrix(rnorm(80L * 40L), nrow = 80L)
  y <- 3 * X[, 1L] - 2 * X[, 20L] + rnorm(80L, sd = .3)
  sampling <- X
  sampling[, seq_len(10L)] <- X[, rep(1L, 10L)]
  missing_y <- y
  missing_y[1L] <- NA_real_
  fixtures <- list(ordinary = list(X = X, y = y),
                   purity_subsampling = list(X = sampling, y = y),
                   invalid_response = list(X = X, y = missing_y))
  for (name in names(fixtures)) {
    fixture <- fixtures[[name]]
    expected <- capture_fit(native, fixture)
    valid <- !identical(name, "invalid_response")
    stopifnot(if (valid) inherits(expected$value, "susie") else
                inherits(expected$value, "susier_test_error"))
    if (identical(name, "purity_subsampling")) {
      stopifnot("The fallback fixture must exercise native RNG use" =
                  !strict_same(expected$rng_before, expected$rng_after))
    }
    callback_calls <- 0L
    callback <- function(fitted, stats) {
      callback_calls <<- callback_calls + 1L
      stopifnot("Native fallback must supply NULL input statistics" = is.null(stats),
                inherits(fitted, "susie"))
      check_namespace()
      fitted
    }
    for (postprocess in list(NULL, callback)) {
      actual <- capture_fit(helpers$runSusieWithFitOptimizations, fixture, postprocess)
      if (!strict_same(expected, actual) ||
          !identical(serialize(expected$value, NULL, version = 3),
                     serialize(actual$value, NULL, version = 3))) {
        stop(name, " fallback must preserve native fit, conditions, bytes, and RNG", call. = FALSE)
      }
      check_namespace()
    }
    stopifnot("Only successful native fits may enter the callback" =
                callback_calls == as.integer(valid))
  }
  message("SuSiE native fallback tests passed: exact fits, RNG, NULL stats, and namespace")
  invisible(TRUE)
}
