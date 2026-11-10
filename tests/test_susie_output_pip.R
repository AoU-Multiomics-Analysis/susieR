# Catch effect-specific alpha values exported as variant PIPs, including the
# loss of native PIPs for variants outside credible sets.
run_susie_output_pip_tests <- function(repo_dir = ".") {
  functions_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH",
                              unset = file.path(repo_dir, "R", "utils"))
  helpers <- new.env(parent = globalenv())
  for (file in c("InitFunctions.R", "SusieFunctions.R", "SusieOutputFunctions.R")) {
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

  # Two correlated pairs produce two credible sets. Keep all fitted effects
  # active so their combined PIPs differ from the assigned effect's alpha.
  set.seed(20261109L)
  n <- 120L
  p <- 20L
  variant_ids <- paste0("chr1_", 1000L + seq_len(p), "_A_G")
  x <- matrix(rnorm(n * p), nrow = n,
              dimnames = list(NULL, variant_ids))
  x[, 2L] <- x[, 1L] + rnorm(n, sd = 0.1)
  x[, 4L] <- x[, 3L] + rnorm(n, sd = 0.1)
  y <- cbind(0.65 * x[, 1L] - 0.55 * x[, 3L] + rnorm(n))
  fitted <- suppressMessages(susieR::susie(
    x, y, L = 10L, estimate_prior_variance = FALSE,
    scaled_prior_variance = 0.1, compute_univariate_zscore = TRUE
  ))
  fitted$variant_id <- variant_ids
  stopifnot("The fixture must produce multiple credible sets" =
              length(fitted$sets$cs) >= 2L)
  cs_variants <- unlist(fitted$sets$cs, use.names = FALSE)
  assigned_effects <- rep(fitted$sets$cs_index, lengths(fitted$sets$cs))
  assigned_alpha <- fitted$alpha[cbind(assigned_effects, cs_variants)]
  stopifnot("The fixture must separate native PIPs from assigned alpha" =
              any(abs(fitted$pip[cs_variants] - assigned_alpha) > 0.01),
            "The fixture must retain variants outside credible sets" =
              length(setdiff(seq_len(p), cs_variants)) > 0L)

  results <- list(gene1 = fitted)
  regions <- tibble::tibble(phenotype_id = "gene1", region = "chr1_1_2000")
  tables <- helpers$formatSusieOutputs(results, regions)
  cached <- helpers$formatSusieOutputs(
    results, regions, extracted_results = list(gene1 = helpers$extractResults(fitted))
  )

  # Use the native fit as the oracle. Do not recalculate PIPs from alpha or
  # credible-set membership in the test.
  failures <- character()
  check <- function(condition, label) {
    if (!isTRUE(condition)) failures <<- c(failures, label)
  }
  same_pips <- function(table) {
    index <- match(table$variant, fitted$variant_id)
    !anyNA(index) && !anyNA(table$pip) &&
      isTRUE(all.equal(unname(table$pip), unname(fitted$pip[index]),
                       tolerance = 0, check.attributes = FALSE))
  }
  check(nrow(tables$in_cs_variant_df) > 0L && same_pips(tables$in_cs_variant_df),
        "Credible-set output PIPs must equal native fitted PIPs, not assigned alpha")
  check(nrow(tables$variant_df) == p &&
          !anyDuplicated(tables$variant_df$variant) &&
          setequal(tables$variant_df$variant, variant_ids),
        "Full output must contain each fitted variant exactly once")
  check(same_pips(tables$variant_df),
        "Full output PIPs must equal native fitted PIPs for every variant")
  outside_cs <- tables$variant_df %>%
    dplyr::filter(is.na(cs_index))
  check(nrow(outside_cs) > 0L && same_pips(outside_cs),
        "Variants outside credible sets must retain native PIPs instead of NA")
  check(identical(tables, cached),
        "Cached extractions must preserve the complete PIP outputs")
  output_dir <- tempfile("susier-native-pip-output-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
  files <- helpers$writeSusieOutputs(results, tables, file.path(output_dir, "gene1"))
  check(same_pips(arrow::read_parquet(files[["summary"]])),
        "Serialized credible-set output must preserve native PIPs")
  check(same_pips(arrow::read_parquet(files[["full"]])),
        "Serialized full output must preserve every native PIP")
  if (length(failures)) stop(paste(failures, collapse = "\n"), call. = FALSE)

  message("SuSiE output PIP tests passed: native PIPs, unassigned variants, and cached output")
  invisible(TRUE)
}
