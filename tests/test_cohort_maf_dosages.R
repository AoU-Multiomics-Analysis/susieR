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
  message("Cohort MAF dosage tests passed: diploid range, folding, sample subset, and predictor order")
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  run_cohort_maf_dosage_tests(if (length(args)) args[[1L]] else ".")
}
