# Frequencies use uncentered diploid ALT dosage, not genotype residuals or SD.
run_cohort_maf_dosage_tests <- function(repo_dir = ".") {
  functions_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH", unset = file.path(repo_dir, "R", "utils"))
  helpers <- new.env(parent = globalenv())
  sys.source(file.path(functions_path, "SusieFunctions.R"), envir = helpers)
  dosage <- rbind(
    c(0, 0, 0, 1), c(2, 2, 2, 1), c(.2, .4, .6, .8), c(0, 0, 0, 0),
    c(2, 2, 2, 2), c(0, 1, 1, 2), c(-1, 1, 1, 1), c(0, 0, 0, 3),
    c(0, 1, NA_real_, 2), c(0, 1, Inf, 2))
  rownames(dosage) <- rep("same_label", nrow(dosage))
  before <- serialize(dosage, NULL)
  expected <- c(.125, .125, .25, 0, 0, .5, rep(NA_real_, 4L))
  actual <- helpers$cohortMafFromDosages(dosage)
  stopifnot("MAF must fold ALT frequencies and preserve physical predictor order" =
              identical(actual, expected),
            "MAF recording must not change the genotype input" =
              identical(serialize(dosage, NULL), before),
            "Empty variant matrices must have a numeric MAF vector" =
              identical(helpers$cohortMafFromDosages(dosage[FALSE, , drop = FALSE]), numeric()))
  # The selected sample set can have a different frequency from the raw cohort.
  subset <- dosage[1:2, 1:3, drop = FALSE]
  stopifnot(identical(helpers$cohortMafFromDosages(subset), c(0, 0)))
  # Counting imputed cells, assuming ALT is minor, or rounding fractional
  # dosages must fail these hand-checked observed-allele fixtures.
  observed <- rbind(
    c(0, NA, -1, 1), c(2, NA, -1, 1), c(.25, NA, 1.5, -1),
    c(NA, -1, NA, -1), c(2, 2, 2, 2), c(0, 0, 0, 0),
    c(-.5, NA, 0, 1), c(3, 0, NA, -1), c(Inf, 1, -1, NA),
    c(NaN, 0, 1, 2))
  rownames(observed) <- rep("same_label", nrow(observed))
  observed_before <- serialize(observed, NULL)
  counts <- helpers$cohortAlleleCountsFromDosages(observed)
  stopifnot("AN must count observed diploid alleles before imputation" =
              identical(counts$cohort_an, c(4, 4, 4, 0, 8, 8, NA_real_, NA_real_, NA_real_, 6)),
            "MAC must fold observed ALT dosage and retain fractional counts" =
              identical(counts$cohort_mac, c(1, 1, 1.75, 0, 0, 0, NA_real_, NA_real_, NA_real_, 3)),
            "Count recording must keep the source matrix unchanged" =
              identical(serialize(observed, NULL), observed_before))
  selected <- helpers$cohortAlleleCountsFromDosages(observed[1:3, c(1, 4), drop = FALSE])
  stopifnot("Counts must use only the selected samples" =
              identical(selected, list(cohort_an = c(4, 4, 2), cohort_mac = c(1, 1, .25))),
            "Empty predictors must have numeric count vectors" =
              identical(helpers$cohortAlleleCountsFromDosages(observed[FALSE, , drop = FALSE]),
                        list(cohort_an = double(), cohort_mac = double())))
  # An unfiltered matrix must retain both physical rows with the same label.
  duplicate_rows <- observed[1:3, , drop = FALSE]
  filtered <- helpers$filterMAF(duplicate_rows, NULL, return_indices = TRUE)
  stopifnot("Filter sidecars must preserve physical rows with repeated labels" =
              identical(filtered, list(genotypes = duplicate_rows, row_indices = 1:3)))
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(stringr)
  })
  sys.source(file.path(functions_path, "ImportFunctions.R"), envir = helpers)
  named_rows <- matrix(0:11, nrow = 3L,
    dimnames = list(c("chrchr1_1_A_G", "chr1_1_A_G", "chr1_2_A_T"), as.character(1:4)))
  allowlist <- tempfile("cohort-count-filter-")
  on.exit(unlink(allowlist), add = TRUE)
  readr::write_tsv(tibble::tibble(variant = c("chr1_2_A_T", "chr1_1_A_G")), allowlist)
  selected_rows <- helpers$filterMAF(named_rows, NULL, variant_list = allowlist,
                                     return_indices = TRUE)
  expected_rows <- named_rows[c(1L, 3L), , drop = FALSE]
  rownames(expected_rows) <- c("chr1_1_A_G", "chr1_2_A_T")
  stopifnot("Allowlist sidecars must follow the legacy first physical row and normalization" =
              identical(selected_rows, list(genotypes = expected_rows, row_indices = c(1L, 3L))))
  ancestry_rows <- rbind(c(0, 0, 2, 2), c(0, 1, 0, 1), c(0, 0, 0, 1))
  dimnames(ancestry_rows) <- list(c("v1", "v2", "v3"), as.character(1:4))
  ancestry <- tibble::tibble(research_id = 1:4,
                            ancestry_pred_other = c("eur", "eur", "afr", "afr"))
  selected_ancestry <- helpers$filterMAF(ancestry_rows, ancestry, MAF_threshold = .2,
                                        return_indices = TRUE)
  stopifnot("Ancestry filter sidecars must retain the selected physical rows" =
              identical(selected_ancestry,
                        list(genotypes = ancestry_rows[2:3, , drop = FALSE], row_indices = 2:3)))
  message("Cohort MAF dosage tests passed: diploid range, folding, sample subset, and predictor order")
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  run_cohort_maf_dosage_tests(if (length(args)) args[[1L]] else ".")
}
