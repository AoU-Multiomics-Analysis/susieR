# Repeated labels must not expand physical predictors or assign a credible set
# to a different predictor with the same label. Native column indices are the
# membership oracle; native fitted values are the numeric oracle.
run_susie_duplicate_variant_label_tests <- function(repo_dir = ".") {
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

  set.seed(20261109L)
  n <- 120L
  p <- 20L
  unique_ids <- paste0("chr1_", 1000L + seq_len(p), "_A_G")
  x <- matrix(rnorm(n * p), nrow = n,
              dimnames = list(NULL, unique_ids))
  x[, 2L] <- x[, 1L] + rnorm(n, sd = 0.1)
  x[, 4L] <- x[, 3L] + rnorm(n, sd = 0.1)
  y <- cbind(0.65 * x[, 1L] - 0.55 * x[, 3L] + rnorm(n))
  regions <- tibble::tibble(phenotype_id = "gene1", region = "chr1_1_2000")
  output_dir <- tempfile("susier-repeated-label-output-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)

  failures <- character()
  check <- function(condition, label) {
    if (!isTRUE(condition)) failures <<- c(failures, label)
  }
  same_values <- function(actual, expected) {
    isTRUE(all.equal(unname(as.matrix(actual)), unname(as.matrix(expected)),
                     tolerance = 0, check.attributes = FALSE))
  }
  posterior_columns <- c(
    "variant_id", "pip", "z", "posterior_mean", "posterior_sd",
    "X_column_scale_factors", "phenotype_residual_SD", paste0("alpha", 1:10),
    paste0("mu_", 1:10), paste0("mu2_", 1:10), "cs_id", "cs_log10bf",
    "cs_avg_r2", "cs_min_r2", "cs_size", "low_purity",
    "has_covariate_explained_variant"
  )
  unique_tables <- NULL
  for (case in c("unique", "both_in_cs", "one_in_cs")) {
    variant_ids <- unique_ids
    if (case == "both_in_cs") variant_ids[[2L]] <- variant_ids[[1L]]
    if (case == "one_in_cs") variant_ids[[20L]] <- variant_ids[[1L]]
    colnames(x) <- variant_ids
    fitted <- suppressMessages(susieR::susie(
      x, y, L = 10L, estimate_prior_variance = FALSE,
      scaled_prior_variance = 0.1, compute_univariate_zscore = TRUE
    ))
    fitted$variant_id <- variant_ids
    stopifnot("The fixture must have two non-overlapping native sets" =
                identical(fitted$sets$cs, list(L1 = c(1L, 2L), L2 = c(3L, 4L))),
              "Repeated labels must refer to distinct genotype columns" =
                !identical(x[, 1L], x[, if (case == "one_in_cs") 20L else 2L]),
              "The repeated predictors must retain distinct fitted values" =
                fitted$pip[[1L]] != fitted$pip[[if (case == "one_in_cs") 20L else 2L]])
    fit_before <- serialize(fitted, NULL)
    seed_before <- .Random.seed
    extracted <- helpers$extractResults(fitted)
    tables <- helpers$formatSusieOutputs(list(gene1 = fitted), regions)
    label <- function(message) paste0(case, ": ", message)

    # The expected memberships use physical column positions, not label matches.
    expected_membership <- c("L1", "L1", "L2", "L2", rep(NA_character_, p - 4L))
    check(nrow(extracted$variant_df) == p,
          label("export must contain one row per physical predictor"))
    check(identical(extracted$variant_df$variant_id, variant_ids),
          label("export must preserve every label and input row order"))
    check(identical(extracted$variant_df$cs_id, expected_membership),
          label("credible-set membership must follow native predictor indices"))
    check(identical(names(extracted$variant_df), posterior_columns),
          label("public extraction columns and their order must stay unchanged"))
    check(same_values(extracted$variant_df$pip, fitted$pip) &&
            same_values(extracted$variant_df$z, fitted$z[, 1L]) &&
            same_values(extracted$variant_df$posterior_mean,
                        susieR::susie_get_posterior_mean(fitted)) &&
            same_values(extracted$variant_df$posterior_sd,
                        susieR::susie_get_posterior_sd(fitted)) &&
            same_values(extracted$variant_df$X_column_scale_factors,
                        fitted$X_column_scale_factors),
          label("each physical row must retain its native PIP, z, mean, SD, and scale"))
    for (matrix_name in c("alpha", "mu", "mu2")) {
      columns <- if (matrix_name == "alpha") paste0("alpha", 1:10) else {
        paste0(matrix_name, "_", 1:10)
      }
      check(same_values(extracted$variant_df[, columns], t(fitted[[matrix_name]])),
            label(paste(matrix_name, "must retain native values in predictor order")))
    }
    check(identical(extracted$lbf_df$variant_id, variant_ids) &&
            same_values(extracted$lbf_df[, -1L], t(fitted$lbf_variable)),
          label("LBF output must retain one native row per physical predictor"))
    check(identical(extracted$cs_df$cs_size, c(2L, 2L)) &&
            same_values(extracted$cs_df$cs_avg_r2,
                        fitted$sets$purity$mean.abs.corr^2) &&
            same_values(extracted$cs_df$cs_min_r2,
                        fitted$sets$purity$min.abs.corr^2),
          label("native credible-set sizes and purity must remain unchanged"))
    check(nrow(tables$variant_df) == p && nrow(tables$in_cs_variant_df) == 4L &&
            identical(tables$in_cs_variant_df$variant, variant_ids[1:4]) &&
            identical(tables$variant_df$cs_index, expected_membership) &&
            same_values(tables$in_cs_variant_df$pip, fitted$pip[1:4]),
          label("full and summary tables must preserve physical rows and assigned values"))

    if (case == "unique") {
      unique_tables <- tables
    } else {
      # The input labels are the only difference between these native fits.
      # Compare all public table columns, including their order and attributes.
      expected <- unique_tables
      expected$variant_df$variant <- variant_ids
      expected$in_cs_variant_df$variant <- variant_ids[1:4]
      expected$lbf_df$variant <- variant_ids
      if (case == "both_in_cs") {
        expected$variant_df$position[[2L]] <- "1001"
        expected$in_cs_variant_df$position[[2L]] <- 1001L
        expected$lbf_df$position[[2L]] <- 1001L
      } else {
        expected$variant_df$position[[20L]] <- "1001"
        expected$lbf_df$position[[20L]] <- 1001L
      }
      for (table in c("variant_df", "in_cs_variant_df")) {
        for (column in intersect(names(expected[[table]]),
                                 c("pip", "posterior_mean", "posterior_sd",
                                   "X_column_scale_factors"))) {
          if (!is.null(names(expected[[table]][[column]]))) {
            names(expected[[table]][[column]]) <- expected[[table]]$variant
          }
        }
      }
      check(identical(tables, expected),
            label("all public outputs must equal the unique-label outputs after label changes"))
    }

    files <- helpers$writeSusieOutputs(list(gene1 = fitted), tables,
                                      file.path(output_dir, case))
    output_tables <- c(summary = "in_cs_variant_df", full = "variant_df", lbf = "lbf_df")
    for (path_key in names(output_tables)) {
      file <- output_tables[[path_key]]
      check(identical(as.data.frame(arrow::read_parquet(files[[path_key]])),
                      as.data.frame(tables[[file]])),
            label(paste(file, "must survive the Parquet write unchanged")))
    }
    check(identical(serialize(fitted, NULL), fit_before) &&
            identical(readRDS(files[["rds"]]), list(gene1 = fitted)),
          label("export must not change the native fit or saved fit"))
    check(identical(.Random.seed, seed_before),
          label("export must not change the RNG state"))
  }

  # Keep native numeric values and set purity, but control CS membership to
  # cover two equal labels in separate sets and a true physical overlap.
  overlapping <- fitted
  overlapping$sets$cs <- list(L1 = c(1L, 2L), L2 = c(20L, 4L), L3 = c(2L, 5L))
  overlapping$sets$cs_index <- 1:3
  overlapping$sets$purity <- fitted$sets$purity[c(1L, 2L, 2L), , drop = FALSE]
  rownames(overlapping$sets$purity) <- c("L1", "L2", "L3")
  overlap_extracted <- helpers$extractResults(overlapping)
  expected_membership <- rep(NA_character_, p)
  expected_membership[c(1L, 2L)] <- "L1"
  expected_membership[c(20L, 4L)] <- "L2"
  check(identical(overlap_extracted$variant_df$cs_id, expected_membership),
        "Equal labels in separate sets must retain separate native-index membership")
  check(identical(overlap_extracted$cs_df, extracted$cs_df),
        "Only a shared physical predictor must suppress the later credible set")

  if (length(failures)) stop(paste(failures, collapse = "\n"), call. = FALSE)
  message("SuSiE repeated-label tests passed: native row values, index membership, schemas, and saved outputs")
  invisible(TRUE)
}
