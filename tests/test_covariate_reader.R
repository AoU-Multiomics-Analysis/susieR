# Compare the production reader with the original parser and matrix conversion.
run_covariate_reader_tests <- function(repo_dir = ".") {
  helper_env <- new.env(parent = globalenv())
  helper_env$`%>%` <- magrittr::`%>%`
  helper_env$pull <- dplyr::pull
  functions_path <- Sys.getenv(
    "SUSIER_FUNCTIONS_PATH", unset = file.path(repo_dir, "R", "utils")
  )
  sys.source(file.path(functions_path, "ImportFunctions.R"), envir = helper_env)
  stopifnot(exists("readQtlCovariateTable", envir = helper_env, inherits = FALSE))

  fixture_dir <- tempfile("covariate-reader-")
  dir.create(fixture_dir)
  on.exit(unlink(fixture_dir, recursive = TRUE), add = TRUE)
  original_reader <- function(path) {
    utils::read.table(path, check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)
  }
  # Retain the complete original matrix conversion as the independent reference.
  original_matrix <- function(path, load_colnames) {
    pc_matrix <- original_reader(path)
    col_names <- dplyr::pull(pc_matrix, ID)
    pc_transpose <- t(pc_matrix[, -1])
    colnames(pc_transpose) <- pc_matrix$SampleID
    pc_df <- dplyr::mutate(as.data.frame(pc_transpose), genotype_id = rownames(pc_transpose)) |>
      dplyr::as_tibble() |>
      dplyr::select(genotype_id, dplyr::everything())
    pc_matrix <- as.matrix(dplyr::select(pc_df, -genotype_id))
    rownames(pc_matrix) <- pc_df$genotype_id
    if (load_colnames) colnames(pc_matrix) <- col_names
    pc_matrix
  }
  capture <- function(expr) {
    warnings <- character()
    value <- tryCatch(
      withCallingHandlers(expr, warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }),
      error = function(condition) condition
    )
    if (inherits(value, "error")) {
      list(error = conditionMessage(value), class = class(value), warnings = warnings)
    } else {
      list(value = value, warnings = warnings)
    }
  }
  fixtures <- list(
    numeric = c("ID\ts2\ts1", "pc2\t0.12345678901234567\t-1.2345678901234567", "pc1\t1e-30\tNaN", "pc3\tInf\t-Inf"),
    integer = c("ID\ts2\ts1", "pc2\t1\t-2", "pc1\t0\t2147483648"),
    logical = c("ID\ts2\ts1", "pc2\tTRUE\tF", "pc1\tFALSE\tT", "pc3\tNA\tNA"),
    all_missing = c("ID\ts2\ts1", "pc2\tNA\tNA", "pc1\tNA\tNA"),
    numeric_ids = c("ID\t002\t001", "02\t01\t2", "01\t2\t3"),
    duplicate_headers = c("ID\ts1\ts1", "pc2\t1\t2", "pc1\t3\t4"),
    literal_na_header = c("ID\tNA\ts2", "pc1\t1\t2", "pc2\t3\t4"),
    quoted = c('"ID"\t"s 2"\t"s1"', '"pc 2"\t1\t2', '"pc 1"\t3\t4'),
    inline_comment = c("ID\ts2\ts1", "pc2\t1\t2 # first", "pc1\t3\t4 # second"),
    whitespace = c("ID s2 s1", "pc2 1 2", "pc1 3 4"),
    blank_line = c("ID\ts2\ts1", "pc2\t1\t2", "", "pc1\t3\t4"),
    empty_field = c("ID\ts2\ts1", "pc2\t\t2", "pc1\t3\t4"),
    missing_field = c("ID\ts2\ts1", "pc2\t1", "pc1\t3\t4"),
    full_10000_samples_20_covariates = c(
      paste(c("ID", paste0("s", seq_len(10000L))), collapse = "\t"),
      purrr::map_chr(seq_len(20L), ~ paste(c(paste0("pc", .x),
        format(seq_len(10000L) / (.x + 0.1), digits = 17, trim = TRUE)), collapse = "\t"))
    ),
    implicit_row_names = c("ID\ts2", "pc2\t1\t2", "pc1\t3\t4"),
    quoted_multiline = c('ID\ts2\ts1', '"pc', '2"\t1\t2', 'pc1\t3\t4'),
    header_only = "ID\ts2\ts1",
    empty = character(),
    utf8 = c("ID\tsé\ts1", "pç2\t1\t2", "pc1\t3\t4"),
    spaces_in_fields = c("ID\ts2\ts1", "pc2\t 1 \t2", "pc1\t3\t4"),
    vertical_whitespace = c("ID\ts2\ts1", "pc2\t1\v2\t3", "pc1\t3\t4")
  )
  paths <- purrr::imap(fixtures, function(lines, name) {
    path <- file.path(fixture_dir, paste0(name, ".tsv"))
    writeLines(lines, path, useBytes = TRUE)
    path
  })
  paths$gzip <- file.path(fixture_dir, "numeric.tsv.gz")
  connection <- gzfile(paths$gzip, "wt")
  writeLines(fixtures$numeric, connection)
  close(connection)
  paths$bom <- file.path(fixture_dir, "bom.tsv")
  writeBin(c(as.raw(c(239, 187, 191)), charToRaw(paste0(paste(fixtures$numeric, collapse = "\n"), "\n"))), paths$bom)
  paths$latin1_bytes <- file.path(fixture_dir, "latin1.tsv")
  writeBin(c(charToRaw("ID\ts"), as.raw(233), charToRaw("\ts1\npc2\t1\t2\npc1\t3\t4\n")), paths$latin1_bytes)
  paths$utf16_bytes <- file.path(fixture_dir, "utf16.tsv")
  writeBin(c(as.raw(c(255, 254)), charToRaw("I"), as.raw(0), charToRaw("D"), as.raw(0)), paths$utf16_bytes)
  paths$no_final_newline <- file.path(fixture_dir, "no-final-newline.tsv")
  writeBin(charToRaw(paste(fixtures$numeric, collapse = "\n")), paths$no_final_newline)
  paths$missing_path <- file.path(fixture_dir, "missing.tsv")

  results <- purrr::imap_dfr(paths, function(path, name) {
    stopifnot(identical(capture(helper_env$readQtlCovariateTable(path)), capture(original_reader(path))))
    purrr::walk(c(FALSE, TRUE), function(load_colnames) {
      stopifnot(identical(
        capture(helper_env$importQtlmapCovariates(path, load_colnames)),
        capture(original_matrix(path, load_colnames))
      ))
    })
    data.frame(case = name, exact_table = TRUE, exact_matrix = TRUE)
  })
  invisible(results)
}
