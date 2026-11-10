# This fails if CV output replaces native PIPs with assigned-effect alpha or
# drops native PIPs for fitted variants outside credible sets.
run_susie_cv_output_pip_tests <- function(repo_dir = ".") {
  functions_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH",
                              unset = file.path(repo_dir, "R", "utils"))
  helpers <- new.env(parent = globalenv())
  helpers$`%>%` <- magrittr::`%>%`
  for (file in c("InitFunctions.R", "SusieFunctions.R", "SusieCVFunctions.R")) {
    sys.source(file.path(functions_path, file), envir = helpers)
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(20261109L)
  n <- 120L
  p <- 20L
  variant_ids <- paste0("chr1_", 1000L + seq_len(p), "_A_G")
  x <- matrix(stats::rnorm(n * p), nrow = n,
              dimnames = list(NULL, variant_ids))
  x[, 2L] <- x[, 1L] + stats::rnorm(n, sd = 0.1)
  x[, 4L] <- x[, 3L] + stats::rnorm(n, sd = 0.1)
  y <- cbind(0.65 * x[, 1L] - 0.55 * x[, 3L] + stats::rnorm(n))
  fitted <- suppressMessages(susieR::susie(
    x, y, L = 10L, estimate_prior_variance = FALSE,
    scaled_prior_variance = 0.1, compute_univariate_zscore = TRUE
  ))
  fitted$variant_id <- variant_ids
  cs_variants <- unlist(fitted$sets$cs, use.names = FALSE)
  cs_effects <- rep(fitted$sets$cs_index, lengths(fitted$sets$cs))
  stopifnot("The fixture must produce multiple credible sets" =
              length(fitted$sets$cs) >= 2L,
            "The fixture must separate native PIPs from assigned alpha" =
              any(abs(fitted$pip[cs_variants] -
                        fitted$alpha[cbind(cs_effects, cs_variants)]) > 0.01))

  extracted <- helpers$extractResults(fitted)
  result <- helpers$CleanSusieData(
    purrr::transpose(list(gene1 = extracted)),
    tibble::tibble(phenotype_id = "gene1", region = "chr1_1_2000")
  )
  index <- match(result$variant, fitted$variant_id)
  failures <- character()
  check <- function(condition, label) {
    if (!isTRUE(condition)) failures <<- c(failures, label)
  }
  check(nrow(result) == p && !anyNA(index) && !anyDuplicated(index),
        "CV output must contain each fitted variant exactly once")
  check(isTRUE(all.equal(unname(result$pip), unname(fitted$pip[index]),
                        tolerance = 0, check.attributes = FALSE)),
        "CV output PIPs must equal native fitted PIPs, not assigned alpha")
  outside_cs <- is.na(result$cs_index)
  check(any(outside_cs) && !anyNA(result$pip[outside_cs]),
        "CV variants outside credible sets must retain native PIPs")
  extracted_index <- match(extracted$variant_df$variant_id, fitted$variant_id)
  check("pip" %in% names(extracted$variant_df) &&
          isTRUE(all.equal(unname(extracted$variant_df$pip),
                           unname(fitted$pip[extracted_index]),
                           tolerance = 0, check.attributes = FALSE)),
        "CV extraction must retain native fitted PIPs")
  if (length(failures)) stop(paste(failures, collapse = "\n"), call. = FALSE)

  message("CV output PIP tests passed: native PIPs and unassigned variants")
  invisible(TRUE)
}
