#!/usr/bin/env Rscript

if (!requireNamespace("lintr", quietly = TRUE)) {
  stop("The lintr package is required. Install it with install.packages('lintr').", call. = FALSE)
}

paths <- c("R/scripts", "R/utils")
paths <- paths[dir.exists(paths)]

if (length(paths) == 0) {
  stop("No R source directories found.", call. = FALSE)
}

linters <- list(
  T_and_F_symbol_linter = lintr::T_and_F_symbol_linter(),
  equals_na_linter = lintr::equals_na_linter(),
  missing_argument_linter = lintr::missing_argument_linter(),
  semicolon_linter = lintr::semicolon_linter(),
  seq_linter = lintr::seq_linter(),
  unreachable_code_linter = lintr::unreachable_code_linter(),
  vector_logic_linter = lintr::vector_logic_linter()
)

lints <- unlist(
  lapply(paths, lintr::lint_dir, linters = linters, parse_settings = FALSE),
  recursive = FALSE
)

if (length(lints) > 0) {
  print(lints)
  quit(status = 1)
}

message("R lint ok")
