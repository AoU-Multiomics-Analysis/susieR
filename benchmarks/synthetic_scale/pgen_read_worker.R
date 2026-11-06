# Read one matched region in a fresh R process. File conversion is not timed.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 6L,
          args[[3L]] %in% c("tabix_read", "pgen_read", "pgen_integer"))
.libPaths(c(Sys.getenv("BENCH_LIBRARY"), .libPaths()))
input_dir <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = TRUE)
mode <- args[[3L]]
repo_dir <- normalizePath(args[[4L]], mustWork = TRUE)
expected_dimensions <- rev(as.integer(args[5:6]))

# Load the same namespaces in both processes. The pipeline runs below also
# measure startup. This reader test measures the complete reader call.
for (package in c("eQTLUtils", "pgenlibr", "Rsamtools", "GenomicRanges", "IRanges",
                  "readr", "dplyr", "digest")) {
  stopifnot(requireNamespace(package, quietly = TRUE))
}
source(file.path(repo_dir, "benchmarks", "synthetic_scale", "pgen_region_reader.R"))
invisible(gc())
elapsed <- system.time({
  if (mode == "tabix_read") {
    genotypes <- eQTLUtils::extractGenotypeMatrixFromDosage(
      chr = "1", start = 1, end = 1999999,
      dosage_file = file.path(input_dir, "dosages.tsv.bgz"))
  } else {
    genotypes <- readPgenRegionHardcalls(
      file.path(input_dir, "synthetic.pgen"),
      file.path(input_dir, "synthetic.pvar"),
      file.path(input_dir, "synthetic.psam"),
      chrom = "1", start = 1, end = 1999999,
      as_double = mode != "pgen_integer")
  }
})[["elapsed"]]
stopifnot(identical(dim(genotypes), expected_dimensions))
readr::write_tsv(tibble::tibble(
  mode = mode, reader_seconds = elapsed,
  samples = ncol(genotypes), variants = nrow(genotypes),
  storage_type = typeof(genotypes), matrix_bytes = as.numeric(object.size(genotypes)),
  pgenlibr_version = as.character(utils::packageVersion("pgenlibr"))),
  file.path(output_dir, "reader_metrics.tsv"))

# Check every cell in small blocks after the clock stops. Coerce integer output
# to the common double format so all three readers have the same check values.
blocks <- split(seq_len(nrow(genotypes)), ceiling(seq_len(nrow(genotypes)) / 100L))
checks <- lapply(seq_along(blocks), function(index) {
  values <- genotypes[blocks[[index]], , drop = FALSE]
  storage.mode(values) <- "double"
  observed <- values[!is.na(values)]
  stopifnot(all(observed %in% c(0, 1, 2)))
  tibble::tibble(block = index, first_row = min(blocks[[index]]),
                 last_row = max(blocks[[index]]), missing = sum(is.na(values)),
                 sum_hardcalls = sum(values, na.rm = TRUE),
                 hash = digest::digest(values, algo = "xxhash64"))
})
readr::write_tsv(dplyr::bind_rows(checks), file.path(output_dir, "matrix_checks.tsv"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
message(mode, " read ", nrow(genotypes), " variants and ", ncol(genotypes),
        " samples in ", elapsed, " s")
