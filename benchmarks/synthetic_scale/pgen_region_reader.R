# Benchmark-only adapter for diploid, biallelic hardcalls.
# ReadIntList() returns samples by variants, with NA for missing hardcalls.
# Its indexes are 1-based; NewPgen() requires increasing sample indexes.
# API source: chrchang/plink-ng, 2.0/pgenlibr/R/RcppExports.R.
# Each call includes metadata loading, file opening, ordering and transpose.
readPgenRegionHardcalls <- function(pgen_path,
                                   pvar_path,
                                   psam_path,
                                   chrom,
                                   start,
                                   end,
                                   sample_ids = NULL,
                                   as_double = TRUE) {
  required_packages <- c("pgenlibr", "readr", "dplyr")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages) > 0L) {
    stop("Missing benchmark packages: ", paste(missing_packages, collapse = ", "),
         call. = FALSE)
  }
  input_paths <- c(pgen = pgen_path, pvar = pvar_path, psam = psam_path)
  if (length(input_paths) != 3L || anyNA(input_paths) ||
      any(!file.exists(input_paths)) || any(file.access(input_paths, 4L) != 0L)) {
    stop("PGEN, PVAR and PSAM inputs must be readable local files", call. = FALSE)
  }
  if (length(chrom) != 1L || is.na(chrom) || !nzchar(as.character(chrom))) {
    stop("chrom must identify one chromosome", call. = FALSE)
  }
  valid_coordinate <- function(value) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value >= 1 && value == floor(value)
  }
  if (!valid_coordinate(start) || !valid_coordinate(end) || start > end) {
    stop("start and end must be positive integer bounds, with start <= end",
         call. = FALSE)
  }
  if (!is.logical(as_double) || length(as_double) != 1L || is.na(as_double)) {
    stop("as_double must be TRUE or FALSE", call. = FALSE)
  }
  normalize_chromosome <- function(value) {
    sub("^(chr)+", "", as.character(value), ignore.case = TRUE)
  }

  read_metadata <- function(path, required_columns, label) {
    # Read PLINK whitespace-separated metadata as strings. A literal NA is an ID.
    connection <- gzfile(path, open = "rt")
    on.exit(close(connection), add = TRUE)
    header <- readLines(connection, n = 1L, warn = FALSE)
    skip <- 0L
    while (length(header) == 1L && startsWith(header, "##")) {
      skip <- skip + 1L
      header <- readLines(connection, n = 1L, warn = FALSE)
    }
    if (length(header) != 1L || !startsWith(header, "#")) {
      stop(label, " requires a header beginning with #", call. = FALSE)
    }
    metadata <- readr::read_table(
      path,
      col_types = readr::cols(.default = readr::col_character()),
      na = character(),
      skip = skip,
      comment = "",
      progress = FALSE,
      show_col_types = FALSE
    )
    names(metadata) <- sub("^#", "", names(metadata))
    if (!all(required_columns %in% names(metadata)) ||
        anyDuplicated(names(metadata)) > 0L || nrow(readr::problems(metadata)) > 0L) {
      stop(label, " has invalid columns or records", call. = FALSE)
    }
    metadata
  }

  pvar_metadata <- read_metadata(pvar_path, c("CHROM", "POS", "ID", "REF", "ALT"),
                                 "PVAR")
  psam_metadata <- read_metadata(psam_path, "IID", "PSAM")
  if (nrow(psam_metadata) == 0L || anyNA(psam_metadata$IID) ||
      any(!nzchar(psam_metadata$IID))) {
    stop("PSAM requires nonempty sample IIDs", call. = FALSE)
  }
  if (anyDuplicated(psam_metadata$IID) > 0L) {
    stop("PSAM has duplicate sample IIDs", call. = FALSE)
  }
  positions <- suppressWarnings(as.numeric(pvar_metadata$POS))
  if (anyNA(positions) || any(!is.finite(positions)) || any(positions < 1) ||
      any(positions != floor(positions))) {
    stop("PVAR has invalid variant positions", call. = FALSE)
  }
  # Assign indexes before selecting the region. They refer to the complete PGEN.
  region <- pvar_metadata |>
    dplyr::mutate(
      .pgen_index = dplyr::row_number(),
      .chromosome = normalize_chromosome(CHROM),
      .position = positions
    ) |>
    dplyr::filter(.chromosome == normalize_chromosome(chrom),
                  .position >= start, .position <= end)
  if (nrow(region) == 0L) {
    stop("No variants in the requested PGEN region", call. = FALSE)
  }
  if (any(grepl(",", region$ALT, fixed = TRUE))) {
    stop("The requested PGEN region has multiallelic variants", call. = FALSE)
  }

  if (is.null(sample_ids)) sample_ids <- psam_metadata$IID
  sample_ids <- as.character(sample_ids)
  if (length(sample_ids) == 0L || anyNA(sample_ids) || any(!nzchar(sample_ids))) {
    stop("sample_ids must contain nonempty sample IDs", call. = FALSE)
  }
  if (anyDuplicated(sample_ids) > 0L) {
    stop("sample_ids has duplicate sample IDs", call. = FALSE)
  }
  requested_indexes <- match(sample_ids, psam_metadata$IID)
  if (anyNA(requested_indexes)) {
    stop("Samples are absent from PSAM: ",
         paste(sample_ids[is.na(requested_indexes)], collapse = ", "), call. = FALSE)
  }
  sorted_indexes <- sort(requested_indexes)

  pvar <- NULL
  pgen <- NULL
  on.exit({
    if (!is.null(pgen)) pgenlibr::ClosePgen(pgen)
    if (!is.null(pvar)) pgenlibr::ClosePvar(pvar)
  }, add = TRUE)
  pvar <- pgenlibr::NewPvar(pvar_path)
  if (pgenlibr::GetVariantCt(pvar) != nrow(pvar_metadata)) {
    stop("PVAR parser variant counts do not match", call. = FALSE)
  }
  pgen <- pgenlibr::NewPgen(pgen_path, pvar = pvar,
                           sample_subset = as.integer(sorted_indexes))
  if (pgenlibr::GetRawSampleCt(pgen) != nrow(psam_metadata)) {
    stop("PSAM and PGEN sample counts do not match", call. = FALSE)
  }
  hardcalls <- pgenlibr::ReadIntList(pgen, as.integer(region$.pgen_index))
  sample_order <- match(requested_indexes, sorted_indexes)
  if (!identical(sample_order, seq_len(nrow(hardcalls)))) {
    hardcalls <- hardcalls[sample_order, , drop = FALSE]
  }
  result <- t(hardcalls)
  dimnames(result) <- list(
    paste0("chr", region$.chromosome, "_", region$POS, "_", region$REF, "_", region$ALT),
    sample_ids
  )
  if (as_double) storage.mode(result) <- "double"
  result
}
