# Compare complete matrix checks and summarize fresh-process reader runs.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[[1L]], mustWork = TRUE)
run_names <- c(paste0("tabix_read_", 1:3), paste0("pgen_read_", 1:3), "pgen_integer")
reference <- readr::read_tsv(file.path(root, run_names[[1L]], "matrix_checks.tsv"),
                             show_col_types = FALSE)
# readr attaches a new parser-problem pointer on each read, including rereads
# of the same file. Compare the actual columns and their types instead.
same_columns <- function(left, right) {
  identical(names(left), names(right)) && identical(dim(left), dim(right)) &&
    all(vapply(names(left), function(column) {
      identical(left[[column]], right[[column]])
    }, logical(1)))
}
checks <- list()
metrics <- lapply(run_names, function(run_name) {
  path <- file.path(root, run_name)
  observed <- readr::read_tsv(file.path(path, "matrix_checks.tsv"),
                             show_col_types = FALSE)
  usage <- readr::read_tsv(file.path(path, "resource_usage.tsv"),
                          show_col_types = FALSE)
  result <- readr::read_tsv(file.path(path, "reader_metrics.tsv"),
                           show_col_types = FALSE)
  checks[[length(checks) + 1L]] <<- tibble::tibble(
    run = run_name, all_block_checks_identical = same_columns(reference, observed),
    exit_zero = nrow(usage) == 1L && usage$exit_code == 0L,
    complete_region = nrow(reference) == 200L && result$variants == 20000L &&
      result$samples == 10000L)
  result |>
    dplyr::mutate(run = run_name, process_seconds = usage$wall_seconds,
                  peak_rss_gib = usage$peak_rss_gib,
                  matrix_gib = matrix_bytes / 1024^3)
}) |>
  dplyr::bind_rows()
checks <- dplyr::bind_rows(checks)
summary <- metrics |>
  dplyr::group_by(mode, storage_type) |>
  dplyr::summarise(runs = dplyr::n(),
                   median_reader_seconds = stats::median(reader_seconds),
                   min_reader_seconds = min(reader_seconds),
                   max_reader_seconds = max(reader_seconds),
                   median_peak_rss_gib = stats::median(peak_rss_gib),
                   matrix_gib = stats::median(matrix_gib), .groups = "drop")
readr::write_tsv(checks, file.path(root, "pgen_read_checks.tsv"))
readr::write_tsv(metrics, file.path(root, "pgen_read_runs.tsv"))
readr::write_tsv(summary, file.path(root, "pgen_read_summary.tsv"))
print(checks, n = Inf)
print(summary, width = Inf)
stopifnot(all(checks$all_block_checks_identical), all(checks$exit_zero),
          all(checks$complete_region))
