#!/usr/bin/env Rscript

# Compare two synthetic benchmark runs. Keep this process separate from the
# measured runs so that loading both sets of outputs does not affect peak RSS.

parse_arguments <- function(args) {
  usage <- paste(
    "Usage: Rscript compare_outputs.R --baseline DIR --candidate DIR --output DIR"
  )
  if (length(args) == 1L && args[[1L]] %in% c("-h", "--help")) {
    cat(usage, "\n", sep = "")
    quit(status = 0L)
  }
  if (length(args) != 6L) {
    stop(usage, call. = FALSE)
  }
  keys <- args[seq.int(1L, length(args), by = 2L)]
  values <- args[seq.int(2L, length(args), by = 2L)]
  if (anyDuplicated(keys) || !setequal(keys, c("--baseline", "--candidate", "--output"))) {
    stop(usage, call. = FALSE)
  }
  stats::setNames(as.list(values), sub("^--", "", keys))
}

options <- parse_arguments(commandArgs(trailingOnly = TRUE))
required_packages <- c("arrow", "dplyr", "jsonlite", "purrr", "readr", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0L) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

if (!dir.exists(options$baseline) || !dir.exists(options$candidate)) {
  stop("The baseline and candidate directories must exist.", call. = FALSE)
}
if (!dir.exists(options$output)) {
  dir.create(options$output, recursive = TRUE)
}
if (!dir.exists(options$output)) {
  stop("Cannot create the output directory.", call. = FALSE)
}
options <- purrr::map(options, normalizePath, mustWork = TRUE)

checks <- list()
add_check <- function(artifact, identical, detail = "", phenotype = NA_character_,
                      field = NA_character_, max_abs_diff = NA_real_) {
  checks[[length(checks) + 1L]] <<- tibble::tibble(
    artifact = artifact,
    phenotype = phenotype,
    field = field,
    identical = identical,
    max_abs_diff = max_abs_diff,
    detail = detail
  )
}

describe_difference <- function(baseline, candidate) {
  difference <- tryCatch(
    all.equal(baseline, candidate, tolerance = 0, check.attributes = TRUE),
    error = function(error) paste("Cannot describe difference:", conditionMessage(error))
  )
  if (isTRUE(difference)) {
    return("Values compare equal; identical() detected a type or attribute difference.")
  }
  substr(paste(utils::head(as.character(difference), 8L), collapse = "; "), 1L, 2000L)
}

# Return a diagnostic only. identical() determines whether a check passes.
numeric_max_difference <- function(baseline, candidate) {
  if (is.numeric(baseline) && is.numeric(candidate)) {
    if (length(baseline) != length(candidate)) {
      return(NA_real_)
    }
    if (length(baseline) == 0L) {
      return(0)
    }
    # Drop dimensions and names for this value-only diagnostic. Use doubles so
    # differences between integer RNG states cannot overflow integer arithmetic.
    baseline <- as.double(baseline)
    candidate <- as.double(candidate)
    if (!identical(is.na(baseline), is.na(candidate)) ||
        !identical(is.nan(baseline), is.nan(candidate))) {
      return(Inf)
    }
    finite <- is.finite(baseline) & is.finite(candidate)
    nonfinite <- !finite & !is.na(baseline)
    if (any(baseline[nonfinite] != candidate[nonfinite])) {
      return(Inf)
    }
    if (any(finite)) {
      return(max(abs(baseline[finite] - candidate[finite])))
    }
    return(0)
  }
  if (is.list(baseline) && is.list(candidate)) {
    count <- min(length(baseline), length(candidate))
    if (count == 0L) {
      return(NA_real_)
    }
    differences <- purrr::map_dbl(seq_len(count), function(index) {
      numeric_max_difference(baseline[[index]], candidate[[index]])
    })
    if (all(is.na(differences))) {
      return(NA_real_)
    }
    return(max(differences, na.rm = TRUE))
  }
  NA_real_
}

read_pair <- function(filename, reader) {
  paths <- list(
    baseline = file.path(options$baseline, filename),
    candidate = file.path(options$candidate, filename)
  )
  missing <- names(paths)[!vapply(paths, file.exists, FUN.VALUE = logical(1))]
  if (length(missing) > 0L) {
    add_check(filename, FALSE, paste("Missing file in:", paste(missing, collapse = ", ")))
    return(NULL)
  }
  values <- purrr::imap(paths, function(path, label) {
    tryCatch(
      list(ok = TRUE, value = reader(path)),
      error = function(error) {
        list(ok = FALSE, message = paste(label, conditionMessage(error), sep = ": "))
      }
    )
  })
  failed <- purrr::keep(values, function(value) !value$ok)
  if (length(failed) > 0L) {
    add_check(filename, FALSE, paste(purrr::map_chr(failed, "message"), collapse = "; "))
    return(NULL)
  }
  purrr::map(values, "value")
}

compare_pair <- function(filename, reader) {
  pair <- read_pair(filename, reader)
  if (is.null(pair)) {
    return(NULL)
  }
  same <- identical(pair$baseline, pair$candidate)
  add_check(
    filename,
    same,
    detail = if (same) "Exact match." else describe_difference(pair$baseline, pair$candidate),
    max_abs_diff = numeric_max_difference(pair$baseline, pair$candidate)
  )
  pair
}

fit_summary <- function(fit) {
  if (!is.list(fit)) {
    return(list(class = class(fit), error = "The phenotype result is not a list."))
  }
  list(
    niter = fit$niter,
    converged = fit$converged,
    final_elbo = if (length(fit$elbo) > 0L) utils::tail(fit$elbo, 1L) else NULL,
    credible_set_count = length(fit$sets$cs),
    variant_count = if (!is.null(dim(fit$alpha))) ncol(fit$alpha) else NULL
  )
}

diagnose_fits <- function(pair) {
  if (!is.list(pair$baseline) || !is.list(pair$candidate)) {
    return(list(error = "The saved results must be lists of phenotype results."))
  }
  count <- max(length(pair$baseline), length(pair$candidate))
  purrr::map(seq_len(count), function(index) {
    baseline_present <- index <= length(pair$baseline)
    candidate_present <- index <= length(pair$candidate)
    baseline <- if (baseline_present) pair$baseline[[index]] else NULL
    candidate <- if (candidate_present) pair$candidate[[index]] else NULL
    baseline_name <- names(pair$baseline)[index]
    candidate_name <- names(pair$candidate)[index]
    labels <- c(baseline_name, candidate_name)
    labels <- labels[!is.na(labels) & nzchar(labels)]
    label <- if (length(labels) > 0L) labels[[1L]] else paste0("index_", index)
    same <- baseline_present && candidate_present && identical(baseline, candidate)
    add_check(
      "synthetic_susie.rds", same,
      phenotype = label, field = "<complete phenotype result>",
      detail = if (same) "Exact match." else describe_difference(baseline, candidate),
      max_abs_diff = numeric_max_difference(baseline, candidate)
    )
    if (!same && is.list(baseline) && is.list(candidate)) {
      if (!identical(attributes(baseline), attributes(candidate))) {
        add_check(
          "synthetic_susie.rds", FALSE, phenotype = label, field = "<attributes>",
          detail = describe_difference(attributes(baseline), attributes(candidate))
        )
      }
      fields <- union(names(baseline), names(candidate))
      purrr::walk(fields, function(field) {
        baseline_field_present <- field %in% names(baseline)
        candidate_field_present <- field %in% names(candidate)
        if (!baseline_field_present || !candidate_field_present ||
            !identical(baseline[[field]], candidate[[field]])) {
          add_check(
            "synthetic_susie.rds", FALSE, phenotype = label, field = field,
            detail = if (!baseline_field_present || !candidate_field_present) {
              "The field is missing from one result."
            } else {
              describe_difference(baseline[[field]], candidate[[field]])
            },
            max_abs_diff = numeric_max_difference(baseline[[field]], candidate[[field]])
          )
        }
      })
    }
    list(
      index = index,
      baseline_name = baseline_name,
      candidate_name = candidate_name,
      identical = same,
      baseline = if (baseline_present) fit_summary(baseline) else NULL,
      candidate = if (candidate_present) fit_summary(candidate) else NULL
    )
  })
}

fit_pair <- compare_pair("synthetic_susie.rds", readRDS)
phenotypes <- if (is.null(fit_pair)) NULL else diagnose_fits(fit_pair)
rm(fit_pair)
invisible(gc())

parquet_files <- c(
  "synthetic.parquet",
  "synthetic.lbf_variable.parquet",
  "synthetic.full_susie.parquet"
)
purrr::walk(parquet_files, function(filename) {
  pair <- compare_pair(filename, function(path) arrow::read_parquet(path, as_data_frame = TRUE))
  rm(pair)
  invisible(gc())
})

# Compare lines rather than readr's source-path and column-spec attributes.
invisible(compare_pair(
  "synthetic.variant_position_summary.tsv",
  function(path) readr::read_lines(path, progress = FALSE)
))
purrr::walk(c("rng_before.rds", "rng_after.rds"), function(filename) {
  invisible(compare_pair(filename, readRDS))
})
invisible(compare_pair("environment.rds", readRDS))
invisible(compare_pair("input_fingerprints.tsv", function(path) {
  readr::read_lines(path, progress = FALSE)
}))

json_safe <- function(value) {
  if (is.null(value)) {
    return(NULL)
  }
  if (inherits(value, c("numeric_version", "package_version", "Date", "POSIXt"))) {
    return(as.character(value))
  }
  if (is.data.frame(value)) {
    return(purrr::map(value, json_safe))
  }
  if (is.list(value)) {
    return(purrr::map(unclass(value), json_safe))
  }
  if (is.atomic(value)) {
    if (is.factor(value)) {
      return(as.character(value))
    }
    if (!is.null(names(value))) {
      return(purrr::map(as.list(value), json_safe))
    }
    return(unclass(value))
  }
  paste(capture.output(print(value)), collapse = "\n")
}

read_optional <- function(directory, filename, reader) {
  path <- file.path(directory, filename)
  if (!file.exists(path)) {
    return(list(present = FALSE))
  }
  tryCatch(
    list(present = TRUE, value = json_safe(reader(path))),
    error = function(error) list(present = TRUE, error = conditionMessage(error))
  )
}

optional_pair <- function(filename, reader) {
  list(
    baseline = read_optional(options$baseline, filename, reader),
    candidate = read_optional(options$candidate, filename, reader)
  )
}

read_metrics <- function(path) {
  # Rebuild the tibble from its columns to drop readr file-path attributes.
  data <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  tibble::as_tibble(purrr::map(data, identity))
}

comparison <- dplyr::bind_rows(checks)
passed <- all(comparison$identical)
summary <- list(
  created_at = format(Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ"),
  baseline = options$baseline,
  candidate = options$candidate,
  all_identical = passed,
  checks = nrow(comparison),
  failed_checks = sum(!comparison$identical),
  comparison = comparison,
  phenotypes = phenotypes,
  environment = optional_pair("environment.rds", readRDS),
  stages = optional_pair("stages.tsv", read_metrics),
  resource_usage = optional_pair("resource_usage.tsv", read_metrics),
  comparison_packages = stats::setNames(
    purrr::map_chr(required_packages, function(package) as.character(utils::packageVersion(package))),
    required_packages
  )
)

readr::write_tsv(comparison, file.path(options$output, "comparison.tsv"), na = "NA")
jsonlite::write_json(
  json_safe(summary), file.path(options$output, "summary.json"),
  auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA
)
if (passed) {
  cat("PASS: required outputs, RNG states, input fingerprints, and environments are identical.\n")
} else {
  cat("FAIL: ", sum(!comparison$identical), " of ", nrow(comparison),
      " checks differ or could not be read.\n", sep = "")
}
cat("Report: ", file.path(options$output, "comparison.tsv"), "\n", sep = "")
quit(save = "no", status = if (passed) 0L else 1L)
