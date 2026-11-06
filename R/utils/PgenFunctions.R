# Reusable reader for diploid, biallelic autosomal PGEN hardcalls.
# The reader owns its handles. Call closePgenReader() when work is complete.
# No imputation, dosage conversion or statistical filtering occurs here.

.pgenNormalizeChromosome <- function(value) {
  sub("^(chr)+", "", as.character(value), ignore.case = TRUE)
}

.pgenReadMetadata <- function(path, required_columns, label) {
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
  header_names <- strsplit(sub("^#", "", trimws(header)), "[[:space:]]+")[[1L]]
  if (anyDuplicated(header_names) || !all(required_columns %in% header_names)) {
    stop(label, " has invalid columns or records", call. = FALSE)
  }
  # Read every field as text. Leading zeros and literal NA are valid IDs.
  metadata <- readr::read_table(
    path, col_types = readr::cols(.default = readr::col_character()),
    na = character(), skip = skip, comment = "", progress = FALSE,
    show_col_types = FALSE
  )
  names(metadata) <- sub("^#", "", names(metadata))
  if (!all(required_columns %in% names(metadata)) ||
      anyDuplicated(names(metadata)) > 0L || nrow(readr::problems(metadata)) > 0L) {
    stop(label, " has invalid columns or records", call. = FALSE)
  }
  metadata
}

# Open metadata once and keep one PGEN handle for all subsequent window reads.
# variants$variant_id is the legacy key; original_id is the literal PVAR ID.
# Only requested autosomal windows are supported. Other records may be present.
openPgenReader <- function(pgen_path, pvar_path, psam_path, sample_ids = NULL) {
  input_paths <- list(pgen_path, pvar_path, psam_path)
  valid_path <- function(path) {
    is.character(path) && length(path) == 1L && !is.na(path) && nzchar(path)
  }
  if (!all(vapply(input_paths, valid_path, logical(1)))) {
    stop("PGEN, PVAR and PSAM inputs must be readable local files", call. = FALSE)
  }
  paths <- unlist(input_paths, use.names = FALSE)
  if (any(grepl("^[A-Za-z][A-Za-z0-9+.-]*://", paths))) {
    stop("PGEN input localization error: unresolved URI; supply local files",
         call. = FALSE)
  }
  if (any(!file.exists(paths)) || any(dir.exists(paths)) ||
      any(file.access(paths, 4L) != 0L)) {
    stop("PGEN, PVAR and PSAM inputs must be readable local files", call. = FALSE)
  }
  if (grepl("\\.zst$", pvar_path, ignore.case = TRUE)) {
    stop("Zstandard PVAR files (.pvar.zst) are not supported; ",
         "supply a plain-text PVAR file", call. = FALSE)
  }
  required_packages <- c("pgenlibr", "readr", "tibble")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages)) {
    stop("Missing PGEN reader packages: ", paste(missing_packages, collapse = ", "),
         call. = FALSE)
  }
  pvar_metadata <- .pgenReadMetadata(
    pvar_path, c("CHROM", "POS", "ID", "REF", "ALT"), "PVAR"
  )
  psam_metadata <- .pgenReadMetadata(psam_path, "IID", "PSAM")
  if (!nrow(psam_metadata) || anyNA(psam_metadata$IID) ||
      any(!nzchar(psam_metadata$IID))) {
    stop("PSAM requires nonempty sample IIDs", call. = FALSE)
  }
  if (anyDuplicated(psam_metadata$IID)) {
    stop("PSAM has duplicate sample IIDs", call. = FALSE)
  }
  positions <- suppressWarnings(as.numeric(pvar_metadata$POS))
  chromosomes <- .pgenNormalizeChromosome(pvar_metadata$CHROM)
  if (!nrow(pvar_metadata) || anyNA(positions) || any(!is.finite(positions)) ||
      any(positions < 1) || any(positions != floor(positions))) {
    stop("PVAR has invalid variant positions", call. = FALSE)
  }
  if (anyNA(chromosomes) || any(!nzchar(chromosomes)) ||
      anyNA(pvar_metadata$REF) || any(!nzchar(pvar_metadata$REF)) ||
      anyNA(pvar_metadata$ALT) || any(!nzchar(pvar_metadata$ALT))) {
    stop("PVAR has invalid chromosomes or alleles", call. = FALSE)
  }
  if (is.null(sample_ids)) sample_ids <- psam_metadata$IID
  if (!is.character(sample_ids) || !length(sample_ids) || anyNA(sample_ids) ||
      any(!nzchar(sample_ids))) {
    stop("sample_ids must contain nonempty sample ID strings", call. = FALSE)
  }
  if (anyDuplicated(sample_ids)) {
    stop("sample_ids has duplicate sample IDs", call. = FALSE)
  }
  requested_indices <- match(sample_ids, psam_metadata$IID)
  if (anyNA(requested_indices)) {
    stop("Samples are absent from PSAM: ",
         paste(sample_ids[is.na(requested_indices)], collapse = ", "), call. = FALSE)
  }
  sorted_indices <- sort(requested_indices)
  reader <- new.env(parent = emptyenv())
  class(reader) <- "PgenReader"
  reader$closed <- FALSE
  reader$pvar <- NULL
  reader$pgen <- NULL
  opened <- FALSE
  on.exit(if (!opened) closePgenReader(reader), add = TRUE)
  reader$pvar <- pgenlibr::NewPvar(pvar_path)
  if (pgenlibr::GetVariantCt(reader$pvar) != nrow(pvar_metadata)) {
    stop("PVAR parser variant counts do not match", call. = FALSE)
  }
  # Check the full sample count before pgenlibr applies the sample subset.
  reader$pgen <- tryCatch(
    pgenlibr::NewPgen(pgen_path, pvar = reader$pvar,
                      raw_sample_ct = nrow(psam_metadata),
                      sample_subset = as.integer(sorted_indices)),
    error = function(error) {
      message <- conditionMessage(error)
      if (grepl("raw_sample_ct", message, fixed = TRUE)) {
        stop("PSAM and PGEN sample counts do not match: ", message, call. = FALSE)
      }
      stop(message, call. = FALSE)
    }
  )
  if (pgenlibr::GetRawSampleCt(reader$pgen) != nrow(psam_metadata)) {
    stop("PSAM and PGEN sample counts do not match", call. = FALSE)
  }
  if (pgenlibr::GetVariantCt(reader$pgen) != nrow(pvar_metadata)) {
    stop("PVAR and PGEN variant counts do not match", call. = FALSE)
  }
  reader$variants <- tibble::tibble(
    pgen_index = seq_len(nrow(pvar_metadata)), chromosome = chromosomes,
    position = positions,
    variant_id = paste0("chr", chromosomes, "_", pvar_metadata$POS, "_",
                        pvar_metadata$REF, "_", pvar_metadata$ALT),
    REF = pvar_metadata$REF, ALT = pvar_metadata$ALT,
    original_id = pvar_metadata$ID
  )
  reader$sample_ids <- sample_ids
  # Keep the full cohort visible even when genotype reads use a sample subset.
  reader$raw_sample_ids <- psam_metadata$IID
  reader$raw_sample_count <- nrow(psam_metadata)
  reader$sample_order <- match(requested_indices, sorted_indices)
  # Sort metadata indices only. Never load or sort a chromosome genotype matrix.
  # Binary searches select positions; global row indices restore PGEN file order.
  reader$position_lookup <- split(reader$variants$pgen_index, chromosomes) |>
    lapply(function(indices) {
      indices <- indices[order(positions[indices], indices)]
      list(positions = positions[indices], indices = indices)
    })
  reg.finalizer(reader, closePgenReader, onexit = TRUE)
  opened <- TRUE
  reader
}

# Bounds are inclusive. The output is always a variants-by-samples matrix.
# missing_code is exactly "NA" or "-1"; integer storage is the default.
readPgenWindow <- function(reader, chrom, start, end,
                           as_double = FALSE, missing_code = "NA") {
  if (!is.environment(reader) || !inherits(reader, "PgenReader")) {
    stop("reader must be a PgenReader", call. = FALSE)
  }
  if (isTRUE(reader$closed)) stop("PgenReader is closed", call. = FALSE)
  if (length(chrom) != 1L || is.na(chrom) ||
      !nzchar(.pgenNormalizeChromosome(chrom))) {
    stop("chrom must identify one chromosome", call. = FALSE)
  }
  chromosome <- .pgenNormalizeChromosome(chrom)
  if (!chromosome %in% as.character(seq_len(22L))) {
    stop("Unsupported ploidy: only diploid autosomal chromosomes 1-22 are supported; ",
         "sex, mitochondrial and other chromosome codes are not supported",
         call. = FALSE)
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
  if (!is.character(missing_code) || length(missing_code) != 1L ||
      is.na(missing_code) || !missing_code %in% c("NA", "-1")) {
    stop("missing_code must be \"NA\" or \"-1\"", call. = FALSE)
  }
  lookup <- reader$position_lookup[[chromosome]]
  if (is.null(lookup)) {
    stop("No variants in the requested PGEN region", call. = FALSE)
  }
  left <- findInterval(start, lookup$positions, left.open = TRUE) + 1L
  right <- findInterval(end, lookup$positions)
  if (left > right) {
    stop("No variants in the requested PGEN region", call. = FALSE)
  }
  indices <- sort(lookup$indices[seq.int(left, right)])
  region <- reader$variants[indices, , drop = FALSE]
  if (any(grepl(",", region$ALT, fixed = TRUE))) {
    stop("The requested PGEN region has multiallelic variants", call. = FALSE)
  }
  hardcalls <- pgenlibr::ReadIntList(reader$pgen, as.integer(indices))
  result <- t(hardcalls[reader$sample_order, , drop = FALSE])
  dimnames(result) <- list(region$variant_id, reader$sample_ids)
  if (missing_code == "-1") result[is.na(result)] <- -1L
  if (as_double) storage.mode(result) <- "double"
  result
}

# Close both owned handles, including on a partial open. Repeated calls are safe.
closePgenReader <- function(reader) {
  if (!is.environment(reader) || !inherits(reader, "PgenReader")) {
    stop("reader must be a PgenReader", call. = FALSE)
  }
  if (isTRUE(reader$closed)) return(invisible(NULL))
  reader$closed <- TRUE
  on.exit({
    reader$pgen <- NULL
    reader$pvar <- NULL
  }, add = TRUE)
  tryCatch({
    if (!is.null(reader$pgen)) pgenlibr::ClosePgen(reader$pgen)
  }, finally = {
    if (!is.null(reader$pvar)) pgenlibr::ClosePvar(reader$pvar)
  })
  invisible(NULL)
}
