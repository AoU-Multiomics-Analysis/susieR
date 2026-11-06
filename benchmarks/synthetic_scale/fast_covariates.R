# Benchmark-only covariate reader. Keep the caller's transpose and dimnames.
# Parse standard TSV fields as text, then use the same base type conversion as
# read.table(). This avoids small changes from fread() numeric conversion.
fast_read_covariates <- function(path) {
  read_original <- function() {
    read.table(
      path,
      check.names = FALSE,
      header = TRUE,
      stringsAsFactors = FALSE
    )
  }

  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  # read.table() accepts whitespace, comments, and quoted fields. Use it for
  # these cases instead of changing their parsing rules in this prototype.
  if (length(lines) == 0L ||
      any(grepl("['\"# ]", lines)) ||
      any(grepl("\\t\\t|^\\t|\\t$", lines))) {
    return(read_original())
  }
  field_counts <- stringr::str_count(lines, stringr::fixed("\t"))
  if (field_counts[[1L]] == 0L || any(field_counts != field_counts[[1L]])) {
    return(read_original())
  }

  parse_warning <- FALSE
  parsed <- tryCatch(
    withCallingHandlers(
      data.table::fread(
        path,
        sep = "\t",
        header = TRUE,
        colClasses = "character",
        na.strings = "NA",
        check.names = FALSE,
        data.table = FALSE
      ),
      warning = function(condition) {
        parse_warning <<- TRUE
        invokeRestart("muffleWarning")
      }
    ),
    error = function(condition) NULL
  )
  if (is.null(parsed) || parse_warning) {
    return(read_original())
  }

  parsed[] <- purrr::map(parsed, ~ type.convert(.x, as.is = TRUE))
  parsed
}
