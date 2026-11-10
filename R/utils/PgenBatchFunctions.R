# Plan complete genes and combine explicitly localized compact batch outputs.
# No helper in this file opens paths stored inside a manifest table.

pgenLocalFile <- function(path, label = "input") {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop(label, " must be one nonempty local path", call. = FALSE)
  }
  if (grepl("^[A-Za-z][A-Za-z0-9+.-]*://", path)) {
    stop(label, " localization error: unresolved URI: ", path, call. = FALSE)
  }
  if (!file.exists(path) || dir.exists(path) || file.access(path, 4L) != 0L) {
    stop(label, " is not a readable local file: ", path, call. = FALSE)
  }
  path
}

pgenCheckIds <- function(values, label, unique = FALSE) {
  if (!is.character(values) || anyNA(values) || any(!nzchar(values)) ||
      any(grepl("[\t\r\n]", values))) {
    stop(label, " must contain nonempty text IDs", call. = FALSE)
  }
  if (unique && anyDuplicated(values)) stop(label, " has duplicate IDs", call. = FALSE)
}

pgenAutosome <- function(chromosome) {
  stripped <- sub("^chr", "", chromosome)
  if (anyNA(stripped) || any(!stripped %in% as.character(1:22))) {
    stop("Only autosome chromosome names 1-22 or chr1-chr22 are supported", call. = FALSE)
  }
  as.integer(stripped)
}

.validatePgenManifest <- function(manifest) {
  columns <- c("chromosome", "pgen", "pvar", "psam")
  if (!is.data.frame(manifest) || !identical(names(manifest), columns)) {
    stop("PGEN manifest schema must be chromosome, pgen, pvar, psam in that order", call. = FALSE)
  }
  if (!nrow(manifest)) stop("PGEN manifest must contain at least one chromosome", call. = FALSE)
  for (column in columns) pgenCheckIds(manifest[[column]], paste("PGEN manifest", column))
  chromosomes <- pgenAutosome(manifest$chromosome)
  if (anyDuplicated(chromosomes)) {
    stop("PGEN manifest contains duplicate chromosome names or aliases", call. = FALSE)
  }
  for (column in c("pgen", "pvar", "psam")) {
    # Inspect path text only. WDL promotes it to File before batch localization.
    values <- manifest[[column]]
    absolute <- startsWith(values, "/") |
      grepl("^[A-Za-z][A-Za-z0-9+.-]*://[^/[:space:]]+(/|$)", values)
    if (any(!absolute)) {
      stop("PGEN manifest ", column, " paths must be absolute cloud or local paths", call. = FALSE)
    }
  }
  tibble::as_tibble(manifest, .name_repair = "minimal")
}

# Read the small manifest only; genotype locations are not opened here.
readPgenManifest <- function(manifest_file) {
  manifest <- suppressWarnings(readr::read_tsv(
    pgenLocalFile(manifest_file, "PGEN manifest"),
    col_types = readr::cols(.default = readr::col_character()),
    na = character(), trim_ws = FALSE, name_repair = "minimal", show_col_types = FALSE
  ))
  if (nrow(readr::problems(manifest))) stop("PGEN manifest TSV parse error", call. = FALSE)
  .validatePgenManifest(manifest)
}

# WDL read_tsv() splits fields without CSV quote decoding. Preserve literal text.
writePgenManifest <- function(manifest, output_file) {
  manifest <- .validatePgenManifest(manifest)
  if (!is.character(output_file) || length(output_file) != 1L || is.na(output_file) ||
      !nzchar(output_file) || grepl("[\t\r\n]", output_file) ||
      grepl("^[A-Za-z][A-Za-z0-9+.-]*://", output_file)) {
    stop("Normalized PGEN manifest output must be a local path", call. = FALSE)
  }
  readr::write_tsv(manifest, output_file, col_names = FALSE, quote = "none", escape = "none")
  invisible(output_file)
}

# Group selected input locations without opening the genotype files.
writePgenQcPlan <- function(manifest, batch_chromosomes, output_dir) {
  manifest <- .validatePgenManifest(manifest)
  pgenCheckIds(batch_chromosomes, "Batch chromosomes")
  if (!length(batch_chromosomes)) stop("Batch chromosomes must be nonempty", call. = FALSE)
  manifest_indices <- match(pgenAutosome(batch_chromosomes), pgenAutosome(manifest$chromosome))
  if (anyNA(manifest_indices)) {
    stop("A selected batch chromosome is absent from the PGEN manifest", call. = FALSE)
  }
  # Keep selected manifest row order. All three paths define the cohort.
  unique_pgens <- manifest[sort(unique(manifest_indices)), ] |>
    dplyr::distinct(pgen, pvar, psam, .keep_all = TRUE)
  path_columns <- c("pgen", "pvar", "psam")
  mapping <- manifest[manifest_indices, path_columns] |>
    dplyr::left_join(unique_pgens |>
      dplyr::mutate(pgen_index = dplyr::row_number() - 1L) |>
      dplyr::select(-chromosome), by = path_columns)
  indices <- as.integer(mapping$pgen_index)
  unique_manifest <- file.path(output_dir, "unique_pgens.tsv")
  batch_pgen_indices <- file.path(output_dir, "batch_pgen_indices.json")
  pgen_batch_indices <- file.path(output_dir, "pgen_batch_indices.json")
  batch_within_pgen_indices <- file.path(output_dir, "batch_within_pgen_indices.json")
  groups <- lapply(seq_len(nrow(unique_pgens)) - 1L, function(source_index) {
    which(indices == source_index) - 1L
  })
  within_indices <- integer(length(indices))
  for (group in groups) within_indices[group + 1L] <- seq_along(group) - 1L
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  writePgenManifest(unique_pgens, unique_manifest)
  jsonlite::write_json(indices, batch_pgen_indices, auto_unbox = FALSE)
  jsonlite::write_json(groups, pgen_batch_indices, auto_unbox = FALSE)
  jsonlite::write_json(within_indices, batch_within_pgen_indices, auto_unbox = FALSE)
  invisible(list(unique_manifest = unique_manifest, batch_pgen_indices = batch_pgen_indices,
                 pgen_batch_indices = pgen_batch_indices,
                 batch_within_pgen_indices = batch_within_pgen_indices,
                 unique_pgens = unique_pgens, indices = indices))
}

# These locations describe outputs. Do not open or localize their contents.
writePgenOutputManifest <- function(chromosomes, pgen_uris, pvar_uris, psam_uris, output_file) {
  values <- list(chromosome = chromosomes, pgen = pgen_uris, pvar = pvar_uris, psam = psam_uris)
  lengths <- vapply(values, length, integer(1))
  if (!lengths[[1L]] || any(lengths != lengths[[1L]])) {
    stop("PGEN output metadata lists must have the same nonzero length", call. = FALSE)
  }
  for (name in names(values)) pgenCheckIds(values[[name]], paste("PGEN output metadata", name))
  manifest <- tibble::as_tibble(values) |>
    dplyr::mutate(chromosome = as.character(pgenAutosome(chromosome))) |>
    dplyr::distinct()
  if (anyDuplicated(manifest$chromosome)) {
    stop("Each output chromosome must have one PGEN/PVAR/PSAM triplet", call. = FALSE)
  }
  manifest <- .validatePgenManifest(manifest)
  # Use the same single-line output-path validation as the normalized writer.
  writePgenManifest(manifest, output_file)
  readr::write_tsv(manifest, output_file, quote = "needed")
  invisible(output_file)
}

pgenRelativePath <- function(path) {
  # read_lines task outputs must remain paths relative to the execution directory.
  absolute <- strsplit(normalizePath(path, mustWork = TRUE), "/", fixed = TRUE)[[1]]
  cwd <- strsplit(normalizePath(getwd(), mustWork = TRUE), "/", fixed = TRUE)[[1]]
  common <- 0L
  while (common < min(length(absolute), length(cwd)) &&
         identical(absolute[[common + 1L]], cwd[[common + 1L]])) common <- common + 1L
  paste(c(rep("..", length(cwd) - common), absolute[seq.int(common + 1L, length(absolute))]),
        collapse = "/")
}

pgenWorklistColumns <- function() c(
  "batch_id", "chromosome", "gene_id", "outcome_id", "phenotype_id", "phenotype_pos",
  "requested_window_start", "requested_window_end", "window_start", "window_end", "region"
)

# Decode only the four metadata fields; never parse or rewrite sample values.
.pgenParseBedMetadata <- function(lines, line_numbers) {
  prefixes <- sub("^((?:[^\t]*\t){3}[^\t]*).*", "\\1", lines, perl = TRUE)
  fields <- stringr::str_split_fixed(prefixes, "\t", 4L)
  trimmed <- trimws(as.vector(fields))
  invalid_quote <- grepl('^"', trimmed) &
    !grepl('^"(?:[^"]|"")*"$', trimmed, perl = TRUE)
  if (any(invalid_quote)) {
    i <- ((which(invalid_quote)[[1L]] - 1L) %% length(lines)) + 1L
    stop("Invalid phenotype BED metadata quotes at line ", line_numbers[[i]],
         "; close quoted fields and escape internal quotes as doubled quotes",
         call. = FALSE)
  }
  metadata_names <- c("chromosome", "start", "end", "phenotype_id")
  metadata <- suppressWarnings(readr::read_tsv(
    I(paste(prefixes, collapse = "\n")),
    col_names = metadata_names,
    col_types = readr::cols(.default = readr::col_character()),
    na = character(), trim_ws = TRUE, quote = "\"", comment = "",
    skip_empty_rows = FALSE, name_repair = "minimal",
    num_threads = 1L, progress = FALSE, show_col_types = FALSE
  ))
  if (!identical(names(metadata), metadata_names) ||
      nrow(metadata) != length(lines) || nrow(readr::problems(metadata))) {
    stop("Invalid phenotype BED metadata near line ", line_numbers[[1L]],
         "; use four tab-delimited fields and balanced quotes on one line",
         call. = FALSE)
  }
  tibble::as_tibble(metadata)
}

# Both planning passes use this reader, including compressed inputs of any suffix.
# Only four columns are retained; the callback receives the original full rows.
.pgenReadPhenotypeBed <- function(path, on_chunk = NULL, expected_header = NULL) {
  pgenLocalFile(path, "Phenotype BED")
  connection <- gzfile(path, "rt")
  on.exit(close(connection), add = TRUE)
  line_number <- 0L
  repeat {
    header <- readLines(connection, n = 1L, warn = FALSE)
    line_number <- line_number + length(header)
    if (!length(header)) stop("Missing phenotype BED header", call. = FALSE)
    candidate <- sub("^\ufeff", "", header)
    stripped <- trimws(candidate)
    if (!nzchar(stripped) || startsWith(stripped, "##") ||
        grepl("^(track|browser)([ \t]|$)", stripped)) next
    if (startsWith(stripped, "#") &&
        !grepl("^#(chr|chrom|chromosome)([ \t]|$)", stripped, ignore.case = TRUE)) next
    column_count <- stringr::str_count(candidate, stringr::fixed("\t")) + 1L
    if (column_count < 4L) {
      stop("Phenotype BED header at line ", line_number,
           " must have at least four tab-delimited columns", call. = FALSE)
    }
    header_metadata <- .pgenParseBedMetadata(candidate, line_number)
    chromosome_header <- sub("^#", "", tolower(header_metadata$chromosome))
    if (!chromosome_header %in% c("chr", "chrom", "chromosome") ||
        !tolower(header_metadata$start) %in% c("start", "chromstart") ||
        !tolower(header_metadata$end) %in% c("end", "chromend")) {
      stop("Invalid phenotype BED header at line ", line_number,
           "; expected a chromosome column followed by start, end, and phenotype ID",
           call. = FALSE)
    }
    break
  }
  if (!is.null(expected_header) && !identical(header, expected_header)) {
    stop("Phenotype BED header changed between metadata reading and row splitting",
         call. = FALSE)
  }
  parts <- list()
  row_count <- 0L
  repeat {
    lines <- readLines(connection, n = 1000L, warn = FALSE)
    if (!length(lines)) break
    line_numbers <- line_number + seq_along(lines)
    line_number <- line_number + length(lines)
    keep <- !grepl("^[ \t]*$", lines)
    lines <- lines[keep]
    line_numbers <- line_numbers[keep]
    if (!length(lines)) next
    widths <- stringr::str_count(lines, stringr::fixed("\t")) + 1L
    invalid <- which(widths != column_count)
    if (length(invalid)) {
      i <- invalid[[1L]]
      stop("Invalid phenotype BED row at line ", line_numbers[[i]], ": expected ",
           column_count, " tab-delimited columns but found ", widths[[i]],
           "; quoted tabs and multiline fields are not supported", call. = FALSE)
    }
    metadata <- .pgenParseBedMetadata(lines, line_numbers)
    if (is.null(on_chunk)) {
      parts[[length(parts) + 1L]] <- metadata
    } else {
      on_chunk(metadata, lines)
    }
    row_count <- row_count + length(lines)
  }
  if (!row_count) stop("Phenotype BED has no data rows", call. = FALSE)
  list(header = header, row_count = row_count,
       metadata = if (is.null(on_chunk)) dplyr::bind_rows(parts) else NULL)
}

readPgenPhenotypeMetadata <- function(path) {
  .pgenReadPhenotypeBed(path)$metadata
}

planPgenBatches <- function(phenotype_bed, qtl_results, chromosomes, cis_distance,
                           output_dir, batch_size = 50L, gene_map = NULL,
                           select_top_phenotype_per_cluster = FALSE,
                           top_phenotype_pvalue_column = "qval") {
  pgenLocalFile(phenotype_bed, "Phenotype BED")
  pgenLocalFile(qtl_results, "QTL results")
  if (!is.null(gene_map)) pgenLocalFile(gene_map, "Gene mapping")
  pgenCheckIds(chromosomes, "Supplied chromosomes", unique = TRUE)
  chromosome_numbers <- pgenAutosome(chromosomes)
  if (anyDuplicated(chromosome_numbers)) stop("Ambiguous supplied chromosome names", call. = FALSE)
  if (length(batch_size) != 1L || !is.numeric(batch_size) || is.na(batch_size) ||
      !is.finite(batch_size) || batch_size < 1 || batch_size != floor(batch_size)) {
    stop("batch_size must be a positive integer", call. = FALSE)
  }
  if (length(cis_distance) != 1L || !is.numeric(cis_distance) || is.na(cis_distance) ||
      !is.finite(cis_distance) || cis_distance < 0 || cis_distance != floor(cis_distance)) {
    stop("cis_distance must be a nonnegative integer", call. = FALSE)
  }
  if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) ||
      !nzchar(output_dir) || grepl("^[A-Za-z][A-Za-z0-9+.-]*://", output_dir)) {
    stop("output_dir must be a local directory", call. = FALSE)
  }
  message("Reading the four phenotype BED metadata columns")
  first_pass <- .pgenReadPhenotypeBed(phenotype_bed)
  metadata <- first_pass$metadata
  pgenCheckIds(metadata$phenotype_id, "Phenotype BED", unique = TRUE)
  if (!exists("importQtlmapPermutedPvalues", mode = "function") ||
      !exists("selectTopPhenotypePerCluster", mode = "function")) {
    stop("Source ImportFunctions.R before planning batches", call. = FALSE)
  }
  # Reuse both the existing importer and the existing cluster ranking rule.
  qtl <- importQtlmapPermutedPvalues(qtl_results) %>%
    dplyr::mutate(phenotype_id = as.character(phenotype_id))
  pgenCheckIds(qtl$phenotype_id, "QTL results", unique = TRUE)
  selected <- qtl %>%
    dplyr::filter(phenotype_id %in% metadata$phenotype_id, p_fdr < .05)
  if (isTRUE(select_top_phenotype_per_cluster)) {
    selected <- selectTopPhenotypePerCluster(selected, top_phenotype_pvalue_column)
  }
  if (!nrow(selected)) stop("No phenotypes pass QTL selection", call. = FALSE)
  selected_metadata <- metadata %>% dplyr::filter(phenotype_id %in% selected$phenotype_id)
  if (is.null(gene_map)) {
    selected_metadata <- selected_metadata %>% dplyr::mutate(gene_id = phenotype_id)
  } else {
    mapping <- readr::read_tsv(gene_map, col_types = readr::cols(.default = readr::col_character()),
                               show_col_types = FALSE, name_repair = "minimal")
    if (!identical(names(mapping), c("phenotype_id", "gene_id"))) {
      stop("Gene mapping schema must be phenotype_id, gene_id", call. = FALSE)
    }
    pgenCheckIds(mapping$phenotype_id, "Gene mapping phenotype IDs", unique = TRUE)
    pgenCheckIds(mapping$gene_id, "Gene mapping gene IDs")
    if (any(!selected_metadata$phenotype_id %in% mapping$phenotype_id)) {
      stop("Missing selected phenotype mapping", call. = FALSE)
    }
    selected_metadata <- dplyr::left_join(selected_metadata, mapping, by = "phenotype_id")
  }
  selected_chromosomes <- pgenAutosome(selected_metadata$chromosome)
  chromosome_match <- match(selected_chromosomes, chromosome_numbers)
  if (anyNA(chromosome_match)) {
    stop("Selected chromosome is absent from supplied chromosome names", call. = FALSE)
  }
  selected_metadata <- selected_metadata %>%
    dplyr::mutate(bed_chromosome = chromosome,
                  chromosome = chromosomes[chromosome_match],
                  chromosome_number = selected_chromosomes,
                  start = suppressWarnings(as.numeric(start)),
                  end = suppressWarnings(as.numeric(end)))
  if (any(!is.finite(selected_metadata$start)) || any(!is.finite(selected_metadata$end)) ||
      any(selected_metadata$start < 0) || any(selected_metadata$end < selected_metadata$start) ||
      any(selected_metadata$start != floor(selected_metadata$start)) ||
      any(selected_metadata$end != floor(selected_metadata$end))) {
    stop("Invalid phenotype BED coordinates", call. = FALSE)
  }
  selected_metadata <- selected_metadata %>%
    dplyr::mutate(phenotype_pos = floor((start + end) / 2))
  genes <- selected_metadata %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(chromosome_count = dplyr::n_distinct(chromosome),
                     chromosome = dplyr::first(chromosome),
                     chromosome_number = dplyr::first(chromosome_number),
                     gene_pos = min(phenotype_pos), .groups = "drop")
  if (any(genes$chromosome_count != 1L)) stop("Gene spans more than one chromosome", call. = FALSE)
  genes <- genes %>% dplyr::arrange(chromosome_number, gene_pos, gene_id) %>%
    dplyr::group_by(chromosome) %>%
    dplyr::mutate(chromosome_batch = (dplyr::row_number() - 1) %/% batch_size) %>%
    dplyr::ungroup()
  batch_keys <- genes %>% dplyr::distinct(chromosome, chromosome_batch) %>%
    dplyr::mutate(batch_id = sprintf("batch_%06d", dplyr::row_number()))
  genes <- dplyr::left_join(genes, batch_keys, by = c("chromosome", "chromosome_batch")) %>%
    dplyr::mutate(gene_order = dplyr::row_number())
  worklist <- selected_metadata %>%
    dplyr::left_join(dplyr::select(genes, gene_id, batch_id, gene_order), by = "gene_id") %>%
    dplyr::arrange(gene_order, phenotype_pos, phenotype_id) %>%
    dplyr::transmute(batch_id, chromosome, gene_id, outcome_id = phenotype_id,
                     phenotype_id, phenotype_pos,
                     requested_window_start = phenotype_pos - cis_distance,
                     requested_window_end = phenotype_pos + cis_distance,
                     window_start = pmax(1, requested_window_start),
                     window_end = requested_window_end,
                     region = paste0(bed_chromosome, ":", requested_window_start, "-", requested_window_end))
  if (any(worklist$window_end < worklist$window_start)) stop("Empty cis window", call. = FALSE)
  dir.create(file.path(output_dir, "batches"), recursive = TRUE, showWarnings = FALSE)
  index <- worklist %>% dplyr::group_by(batch_id, chromosome) %>%
    dplyr::summarise(gene_count = dplyr::n_distinct(gene_id), outcome_count = dplyr::n(),
                     .groups = "drop") %>% dplyr::arrange(batch_id)
  worklist_paths <- file.path(output_dir, "batches", paste0(index$batch_id, ".tsv"))
  bed_paths <- file.path(output_dir, "batches", paste0(index$batch_id, ".bed"))
  for (i in seq_len(nrow(index))) {
    readr::write_tsv(dplyr::filter(worklist, batch_id == index$batch_id[[i]]),
                     worklist_paths[[i]], quote = "needed")
  }
  message("Splitting selected phenotype BED rows with the same metadata parser")
  for (path in bed_paths) writeLines(first_pass$header, path, useBytes = TRUE)
  target_batch <- stats::setNames(match(worklist$batch_id, index$batch_id), worklist$phenotype_id)
  written <- stats::setNames(integer(nrow(worklist)), worklist$phenotype_id)
  metadata_offset <- 0L
  split_chunk <- function(chunk, lines) {
    rows <- metadata_offset + seq_len(nrow(chunk))
    if (max(rows) > nrow(metadata) ||
        !identical(chunk, metadata[rows, , drop = FALSE])) {
      stop("Phenotype BED metadata changed between planning and splitting near data row ",
           metadata_offset + 1L, "; check row order, coordinates, and phenotype IDs",
           call. = FALSE)
    }
    metadata_offset <<- metadata_offset + nrow(chunk)
    row_ids <- chunk$phenotype_id
    selected_rows <- which(row_ids %in% names(target_batch))
    if (!length(selected_rows)) return(invisible(NULL))
    groups <- split(selected_rows, target_batch[row_ids[selected_rows]])
    for (batch in names(groups)) {
      destination <- file(bed_paths[[as.integer(batch)]], "at")
      tryCatch(writeLines(lines[groups[[batch]]], destination, useBytes = TRUE),
               finally = close(destination))
    }
    written <<- written + tabulate(match(row_ids[selected_rows], names(written)),
                                   nbins = length(written))
    invisible(NULL)
  }
  second_pass <- .pgenReadPhenotypeBed(phenotype_bed, on_chunk = split_chunk,
                                      expected_header = first_pass$header)
  if (second_pass$row_count != first_pass$row_count) {
    stop("Phenotype BED row count changed between planning and splitting: expected ",
         first_pass$row_count, " rows but read ", second_pass$row_count, call. = FALSE)
  }
  if (any(written != 1L)) {
    stop("Phenotype BED splitting failed: ", sum(written == 0L), " selected IDs missing; ",
         sum(written > 1L), " selected IDs copied more than once. Check the input IDs",
         call. = FALSE)
  }
  index <- index %>% dplyr::mutate(
    worklist = vapply(worklist_paths, pgenRelativePath, character(1)),
    phenotype_bed = vapply(bed_paths, pgenRelativePath, character(1)))
  expected_path <- file.path(output_dir, "expected_outcomes.tsv")
  index_path <- file.path(output_dir, "batch_index.tsv")
  readr::write_tsv(worklist, expected_path, quote = "needed")
  readr::write_tsv(index, index_path, quote = "needed")
  writeLines(index$worklist, file.path(output_dir, "batch_worklists.txt"))
  writeLines(index$phenotype_bed, file.path(output_dir, "batch_phenotype_beds.txt"))
  writeLines(index$chromosome, file.path(output_dir, "batch_chromosomes.txt"))
  # WDL 1.0 has no keys(Map). This output metadata links each batch to the
  # corresponding typed chromosome struct without serializing its File inputs.
  jsonlite::write_json(as.integer(match(index$chromosome, chromosomes) - 1L),
                       file.path(output_dir, "batch_chromosome_indices.json"),
                       auto_unbox = FALSE)
  message("Planned ", nrow(genes), " genes and ", nrow(worklist), " outcomes in ", nrow(index), " batches")
  invisible(list(expected_outcomes = expected_path, batch_index = index_path,
                 worklists = index$worklist, phenotype_beds = index$phenotype_bed,
                 chromosomes = index$chromosome))
}

writePgenBatchRegions <- function(worklist, output_file) {
  if (is.character(worklist)) {
    worklist <- readr::read_tsv(pgenLocalFile(worklist, "Worklist"), show_col_types = FALSE)
  }
  needed <- c("chromosome", "window_start", "window_end")
  if (!all(needed %in% names(worklist)) || !nrow(worklist)) {
    stop("Worklist needs nonempty chromosome/window_start/window_end columns", call. = FALSE)
  }
  pgenAutosome(worklist$chromosome)
  if (!is.numeric(worklist$window_start) || !is.numeric(worklist$window_end) ||
      any(!is.finite(worklist$window_start)) || any(!is.finite(worklist$window_end)) ||
      any(worklist$window_start < 1) || any(worklist$window_end < worklist$window_start) ||
      any(worklist$window_start != floor(worklist$window_start)) ||
      any(worklist$window_end != floor(worklist$window_end))) {
    stop("Invalid inclusive worklist windows", call. = FALSE)
  }
  worklist %>% dplyr::transmute(chromosome, window_start, window_end,
                                label = sprintf("region_%09d", dplyr::row_number())) %>%
    readr::write_tsv(output_file, col_names = FALSE)
  invisible(output_file)
}

pgenReadTextTable <- function(path, numeric_columns = character(), logical_columns = character()) {
  table <- suppressWarnings(readr::read_tsv(pgenLocalFile(path),
                           col_types = readr::cols(.default = readr::col_character()),
                           show_col_types = FALSE, name_repair = "minimal", na = "NA"))
  if (nrow(readr::problems(table))) stop("TSV schema parse error: ", path, call. = FALSE)
  if (anyDuplicated(names(table))) stop("Duplicate schema columns: ", path, call. = FALSE)
  for (column in intersect(numeric_columns, names(table))) {
    value <- suppressWarnings(as.numeric(table[[column]]))
    if (any(is.na(value) & !is.na(table[[column]]))) stop("Invalid numeric column: ", column, call. = FALSE)
    table[[column]] <- value
  }
  for (column in intersect(logical_columns, names(table))) {
    value <- table[[column]]
    if (any(!is.na(value) & !value %in% c("TRUE", "FALSE"))) {
      stop("Invalid boolean column: ", column, call. = FALSE)
    }
    table[[column]] <- as.logical(value)
  }
  table
}

pgenRequireColumns <- function(table, required, label) {
  if (!is.data.frame(table) || anyDuplicated(names(table)) || !all(required %in% names(table))) {
    stop(label, " schema is missing required columns: ", paste(setdiff(required, names(table)), collapse = ", "),
         call. = FALSE)
  }
}

pgenOwnershipKey <- function(table) {
  if (!nrow(table)) return(character())
  # Use length-delimited IDs, so literal separators inside IDs cannot collide.
  columns <- c("batch_id", "gene_id", "outcome_id", "phenotype_id")
  do.call(paste0, lapply(columns, function(column) paste0(nchar(table[[column]], type = "bytes"), ":", table[[column]])))
}

pgenCheckOwnership <- function(table, expected, label, variant_chromosome = FALSE) {
  columns <- c("batch_id", "gene_id", "outcome_id", "phenotype_id")
  pgenRequireColumns(table, columns, label)
  for (column in columns) pgenCheckIds(table[[column]], paste(label, column))
  keys <- pgenOwnershipKey(table)
  match_index <- match(keys, pgenOwnershipKey(expected))
  if (anyNA(match_index)) stop(label, " contains foreign outcome ownership", call. = FALSE)
  if ("molecular_trait_id" %in% names(table) &&
      (anyNA(table$molecular_trait_id) || any(table$molecular_trait_id != table$outcome_id))) {
    stop(label, " molecular_trait_id does not match outcome ownership", call. = FALSE)
  }
  if ("chromosome" %in% names(table) && nrow(table)) {
    chromosome <- expected$chromosome[match_index]
    table_chromosome <- as.character(table$chromosome)
    if (variant_chromosome) {
      # Summary chromosomes describe variants, including gene-owned additional features.
      pgenCheckIds(table$variant, paste(label, "variant ID"))
      if (any(!grepl("^[^_]+_[0-9]+_[^_]+_.+$", table$variant))) {
        stop(label, " has malformed variant IDs", call. = FALSE)
      }
      chromosome <- sub("^(?i:chr)+", "", sub("_.*$", "", table$variant), perl = TRUE)
      table_chromosome <- sub("^(?i:chr)+", "", table_chromosome, perl = TRUE)
    }
    if (any(pgenAutosome(table_chromosome) != pgenAutosome(chromosome))) {
      stop(label, " chromosome does not match expected ",
           if (variant_chromosome) "variant ID" else "outcome", call. = FALSE)
    }
  }
  keys
}

pgenBindCompact <- function(tables, label, empty_compatible = character()) {
  reference <- names(tables[[1L]])
  if (any(vapply(tables, function(table) !setequal(names(table), reference), logical(1)))) {
    stop(label, " has incompatible schema columns across batches", call. = FALSE)
  }
  # Empty legacy tables have numeric IDs. Only zero-row columns can change type.
  nonempty <- which(vapply(tables, nrow, integer(1)) > 0L)
  prototype <- tables[[if (length(nonempty)) nonempty[[1L]] else 1L]]
  for (i in seq_along(tables)) {
    tables[[i]] <- tables[[i]][reference]
    if (!nrow(tables[[i]])) {
      for (column in intersect(empty_compatible, reference)) {
        tables[[i]][[column]] <- prototype[[column]][FALSE]
      }
    }
  }
  family <- function(value) if (is.numeric(value)) "numeric" else typeof(value)
  for (column in reference) {
    types <- vapply(tables, function(table) family(table[[column]]), character(1))
    if (length(unique(types)) != 1L) stop(label, " has incompatible schema type for ", column, call. = FALSE)
  }
  dplyr::bind_rows(tables)
}

# This metadata is optional only for legacy files. Never derive a missing scale.
# A shared state records one value per outcome when full output spans batches.
pgenCheckPhenotypeResidualSD <- function(table, label, seen = NULL) {
  column <- "phenotype_residual_SD"
  if (!column %in% names(table)) table[[column]] <- rep(NA_real_, nrow(table))
  value <- table[[column]]
  if (!is.numeric(value)) {
    stop(label, " requires numeric phenotype_residual_SD", call. = FALSE)
  }
  if (any(!is.na(value) & (!is.finite(value) | value < 0))) {
    stop(label, " phenotype_residual_SD must be nonnegative and finite or missing", call. = FALSE)
  }
  table[[column]] <- as.double(value)
  table[[column]][is.na(table[[column]])] <- NA_real_
  values <- table |>
    dplyr::select(molecular_trait_id, phenotype_residual_SD) |>
    dplyr::distinct()
  if (anyDuplicated(values$molecular_trait_id)) {
    stop(label, " has inconsistent phenotype_residual_SD within an outcome", call. = FALSE)
  }
  if (!is.null(seen)) {
    for (i in seq_len(nrow(values))) {
      outcome <- values$molecular_trait_id[[i]]
      scale <- values$phenotype_residual_SD[[i]]
      if (exists(outcome, envir = seen, inherits = FALSE) &&
          !identical(get(outcome, envir = seen, inherits = FALSE), scale)) {
        stop(label, " has inconsistent phenotype_residual_SD within an outcome", call. = FALSE)
      }
      assign(outcome, scale, envir = seen)
    }
  }
  table
}

# The validated outcome index is the source for old files and unknown fit flags.
pgenCheckFeatureConverged <- function(table, converged, label) {
  column <- "feature_converged"
  if (!column %in% names(table)) table[[column]] <- rep(NA, nrow(table))
  value <- table[[column]]
  if (!is.logical(value)) {
    stop(label, " requires logical feature_converged", call. = FALSE)
  }
  if (any(!is.na(value) & value != converged)) {
    stop(label, " feature_converged does not match the outcome index", call. = FALSE)
  }
  table[[column]] <- dplyr::coalesce(value, converged)
  table
}

aggregatePgenBatches <- function(expected_outcomes, summary_files, credible_set_files,
                                index_files, diagnostic_files, output_dir) {
  message("Checking localized batch files and expected outcomes")
  expected <- pgenReadTextTable(expected_outcomes)
  pgenRequireColumns(expected, pgenWorklistColumns(), "Expected outcomes")
  if (!nrow(expected)) stop("Expected outcomes must not be empty", call. = FALSE)
  for (column in c("batch_id", "chromosome", "gene_id", "outcome_id", "phenotype_id")) {
    pgenCheckIds(expected[[column]], paste("Expected", column))
  }
  pgenAutosome(expected$chromosome)
  expected_keys <- pgenOwnershipKey(expected)
  if (anyDuplicated(expected_keys) || anyDuplicated(expected$outcome_id) || anyDuplicated(expected$phenotype_id)) {
    stop("Expected outcomes contain duplicate IDs", call. = FALSE)
  }
  batch_ids <- unique(expected$batch_id)
  file_sets <- list(summary_files, credible_set_files, index_files, diagnostic_files)
  for (files in file_sets) {
    if (!is.character(files) || length(files) != length(batch_ids) || anyDuplicated(files)) {
      stop("Supply exactly one unique file per expected batch in each file set", call. = FALSE)
    }
    for (path in files) pgenLocalFile(path)
    if (anyDuplicated(normalizePath(files, mustWork = TRUE))) stop("Duplicate localized files", call. = FALSE)
  }
  owner_columns <- c("batch_id", "gene_id", "outcome_id", "phenotype_id")
  index_columns <- c(owner_columns, "chromosome", "gene_key", "seed", "status", "converged", "credible_set_count", "reason")
  summary_columns <- c(owner_columns, "molecular_trait_id", "variant", "chromosome", "position", "ref", "alt",
                       "cs_id", "cs_index", "region", "pip", "z", "cs_min_r2", "cs_avg_r2", "cs_size",
                       "posterior_mean", "posterior_sd", "cs_log10bf", "has_covariate_explained_variant")
  cs_columns <- c(owner_columns, "molecular_trait_id", "cs_id", "cs_index", "region", "cs_log10bf", "cs_avg_r2",
                  "cs_min_r2", "cs_size", "low_purity", "has_covariate_explained_variant")
  diagnostic_columns <- c(owner_columns, "molecular_trait_id", "variant", "input_scale", "residual_sd", "residual_sd_ratio")
  indices <- lapply(index_files, pgenReadTextTable,
                    numeric_columns = c("seed", "credible_set_count", "data_loading_seconds",
                      "data_preparation_seconds", "fine_mapping_seconds", "outcome_total_seconds"),
                    logical_columns = "converged")
  seen_batches <- character()
  for (i in seq_along(indices)) {
    table <- indices[[i]]
    pgenRequireColumns(table, index_columns, "Outcome index")
    if (!nrow(table) || dplyr::n_distinct(table$batch_id) != 1L) {
      stop("Each outcome index must contain one nonempty batch", call. = FALSE)
    }
    seen_batches <- c(seen_batches, table$batch_id[[1L]])
    keys <- pgenCheckOwnership(table, expected, "Outcome index")
    if (anyDuplicated(keys)) stop("Outcome index contains duplicate outcomes", call. = FALSE)
    if (anyNA(table$status) || any(table$status == "FAILED")) stop("FAILED or missing outcome status", call. = FALSE)
    if (any(!table$status %in% c("CONVERGED", "NONCONVERGED")) || anyNA(table$converged) ||
        any(table$converged != (table$status == "CONVERGED"))) {
      stop("Outcome status and converged boolean disagree", call. = FALSE)
    }
    for (column in c("seed", "credible_set_count")) {
      if (any(!is.finite(table[[column]])) || any(table[[column]] < 0) ||
          any(table[[column]] != floor(table[[column]])) ||
          any(table[[column]] > .Machine$integer.max)) stop("Invalid outcome index ", column, call. = FALSE)
    }
    timing_columns <- intersect(c("data_loading_seconds", "data_preparation_seconds",
                                  "fine_mapping_seconds", "outcome_total_seconds"), names(table))
    for (column in timing_columns) {
      value <- table[[column]]
      if (any(!is.na(value) & (!is.finite(value) | value < 0))) {
        stop("Invalid outcome index timing: ", column, call. = FALSE)
      }
    }
    pgenCheckIds(table$gene_key, "Outcome index gene_key")
  }
  if (anyDuplicated(seen_batches) || !setequal(seen_batches, batch_ids)) {
    stop("Duplicate or missing batch file sets", call. = FALSE)
  }
  index <- pgenBindCompact(indices, "Outcome index")
  keys <- pgenOwnershipKey(index)
  if (anyDuplicated(keys)) stop("Duplicate expected outcome statuses", call. = FALSE)
  if (!setequal(keys, expected_keys)) stop("Missing expected outcome statuses", call. = FALSE)
  index <- index[match(expected_keys, keys), ]
  index <- index %>% dplyr::mutate(seed = as.integer(seed), credible_set_count = as.integer(credible_set_count))
  read_parquet <- function(path) tibble::as_tibble(arrow::read_parquet(path))
  summaries <- lapply(summary_files, read_parquet)
  sets <- lapply(credible_set_files, read_parquet)
  diagnostics <- lapply(diagnostic_files, pgenReadTextTable,
                        numeric_columns = c("input_scale", "residual_sd", "residual_sd_ratio"))
  for (i in seq_along(indices)) {
    own_expected <- expected %>% dplyr::filter(batch_id == seen_batches[[i]])
    for (kind in c("summary", "credible sets", "diagnostics")) {
      table <- switch(kind, summary = summaries[[i]], "credible sets" = sets[[i]], diagnostics = diagnostics[[i]])
      required <- switch(kind, summary = summary_columns, "credible sets" = cs_columns, diagnostics = diagnostic_columns)
      pgenRequireColumns(table, required, kind)
      # The only accepted old-schema type differences are on zero-row columns.
      if (!nrow(table)) {
        for (column in intersect(c("molecular_trait_id", "cs_id"), names(table))) table[[column]] <- character()
        for (column in intersect(c("low_purity", "has_covariate_explained_variant"), names(table))) table[[column]] <- logical()
      }
      table_keys <- pgenCheckOwnership(table, own_expected, kind, variant_chromosome = kind == "summary")
      if (kind == "summary") table <- pgenCheckPhenotypeResidualSD(table, "Summary")
      if (kind != "diagnostics") {
        owner <- match(table_keys, pgenOwnershipKey(indices[[i]]))
        table <- pgenCheckFeatureConverged(table, indices[[i]]$converged[owner], kind)
      }
      text_columns <- intersect(c("molecular_trait_id", "variant", "ref", "alt", "cs_id", "region"), names(table))
      for (column in text_columns) pgenCheckIds(table[[column]], paste(kind, "text schema", column))
      if (nrow(table) && "variant" %in% names(table)) {
        row_key <- paste0(pgenOwnershipKey(table), ":", nchar(table$variant), ":", table$variant)
        if (kind == "summary") row_key <- paste0(row_key, ":", nchar(table$cs_id), ":", table$cs_id)
        if (anyDuplicated(row_key)) stop(kind, " contains duplicate rows", call. = FALSE)
      }
      for (column in intersect(c("low_purity", "has_covariate_explained_variant"), names(table))) {
        if (!is.logical(table[[column]]) || anyNA(table[[column]])) {
          stop(kind, " requires a nonmissing logical boolean column: ", column, call. = FALSE)
        }
      }
      numbers <- switch(kind,
                         summary = c("position", "pip", "z", "cs_min_r2", "cs_avg_r2", "cs_size", "posterior_mean", "posterior_sd", "cs_log10bf"),
                         "credible sets" = c("cs_log10bf", "cs_avg_r2", "cs_min_r2", "cs_size"),
                         diagnostics = c("input_scale", "residual_sd", "residual_sd_ratio"))
      if (any(!vapply(table[numbers], is.numeric, logical(1)))) stop(kind, " requires numeric schema columns", call. = FALSE)
      switch(kind, summary = { summaries[[i]] <- table },
             "credible sets" = { sets[[i]] <- table }, diagnostics = { diagnostics[[i]] <- table })
    }
  }
  empty_compatible <- c("molecular_trait_id", "cs_id", "cs_index", "low_purity", "has_covariate_explained_variant")
  summary <- pgenBindCompact(summaries, "Summary", empty_compatible)
  cs <- pgenBindCompact(sets, "Credible sets", empty_compatible)
  diagnostic <- pgenBindCompact(diagnostics, "Diagnostics")
  # Cross-check the count even for valid empty summaries: no-CS must be explicit.
  cs_keys <- pgenOwnershipKey(cs)
  actual_counts <- vapply(expected_keys, function(key) sum(cs_keys == key), integer(1))
  if (any(actual_counts != index$credible_set_count)) {
    stop("Credible-set counts do not match the complete outcome index", call. = FALSE)
  }
  if (nrow(cs) && anyDuplicated(paste0(cs_keys, ":", nchar(cs$cs_id), ":", cs$cs_id))) {
    stop("Duplicate credible set rows", call. = FALSE)
  }
  # The formatter emits every member of each retained high-purity set. Overlapped
  # sets were already removed from both tables by extractResults().
  set_key <- function(table) {
    if (!nrow(table)) return(character())
    paste0(pgenOwnershipKey(table), ":", nchar(table$cs_id, type = "bytes"), ":", table$cs_id)
  }
  cs_set_keys <- set_key(cs)
  summary_set_keys <- set_key(summary)
  if (any(!is.finite(cs$cs_size)) || any(cs$cs_size < 1) || any(cs$cs_size != floor(cs$cs_size))) {
    stop("Credible sets require positive integer cs_size values", call. = FALSE)
  }
  set_match <- match(summary_set_keys, cs_set_keys)
  if (anyNA(set_match)) stop("Summary references an unknown credible set or wrong set ownership", call. = FALSE)
  if (any(cs$low_purity[set_match])) stop("Summary references a low-purity credible set", call. = FALSE)
  if (any(!is.finite(summary$cs_size)) || any(summary$cs_size != cs$cs_size[set_match])) {
    stop("Summary cs_size does not match its credible set", call. = FALSE)
  }
  if (any(summary$has_covariate_explained_variant != cs$has_covariate_explained_variant[set_match])) {
    stop("Summary covariate-explained flag does not match its credible set", call. = FALSE)
  }
  members <- tibble::tibble(set_row = set_match, variant = summary$variant) %>% dplyr::distinct()
  member_counts <- tabulate(members$set_row, nbins = nrow(cs))
  if (any(member_counts[!cs$low_purity] != cs$cs_size[!cs$low_purity])) {
    stop("Summary has missing or excess distinct variants for a high-purity credible set", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  files <- c(summary = file.path(output_dir, "susie.parquet"),
             credible_sets = file.path(output_dir, "credible_sets.parquet"),
             index = file.path(output_dir, "outcome_index.tsv"),
             diagnostics = file.path(output_dir, "covariate_explained_variants.tsv"))
  arrow::write_parquet(summary, files[["summary"]])
  arrow::write_parquet(cs, files[["credible_sets"]])
  readr::write_tsv(index, files[["index"]], quote = "needed")
  readr::write_tsv(diagnostic, files[["diagnostics"]], quote = "needed")
  message("Aggregated ", nrow(index), " outcomes; NONCONVERGED: ", sum(index$status == "NONCONVERGED"))
  invisible(files)
}

# Artifact locations are output metadata. Only the combined outcome index is read.
writeChromosomeGeneManifests <- function(outcome_index_file, gene_rds_uris,
                                        gene_parquet_uris, gene_tsv_uris,
                                        gene_log_uris, gene_provenance_uris,
                                        output_dir) {
  index <- pgenReadTextTable(outcome_index_file, logical_columns = "converged")
  identity_columns <- c("chromosome", "gene_id", "gene_key", "batch_id", "outcome_id", "phenotype_id")
  pgenRequireColumns(index, c(identity_columns, "status", "converged"), "Combined outcome index")
  if (!nrow(index)) stop("Combined outcome index must not be empty", call. = FALSE)
  for (column in identity_columns) pgenCheckIds(index[[column]], paste("Outcome index", column))
  if (anyDuplicated(index$outcome_id) || anyDuplicated(index$phenotype_id)) {
    stop("Combined outcome index contains duplicate outcome or phenotype IDs", call. = FALSE)
  }
  if (anyNA(index$status) || any(!index$status %in% c("CONVERGED", "NONCONVERGED")) ||
      anyNA(index$converged) || any(index$converged != (index$status == "CONVERGED"))) {
    stop("Combined outcome index requires successful statuses and matching converged booleans", call. = FALSE)
  }
  expected_gene_keys <- vapply(index$gene_id, function(gene_id) {
    paste0("gene-", digest::digest(charToRaw(enc2utf8(gene_id)), algo = "sha256", serialize = FALSE))
  }, character(1), USE.NAMES = FALSE)
  if (any(index$gene_key != expected_gene_keys)) {
    stop("Outcome index gene_key does not match the original gene_id", call. = FALSE)
  }
  genes <- index %>%
    dplyr::mutate(chromosome = as.character(pgenAutosome(chromosome))) %>%
    dplyr::distinct(chromosome, gene_id, gene_key, batch_id)
  if (anyDuplicated(genes$gene_id)) {
    stop("Each gene must have exactly one chromosome and batch in the outcome index", call. = FALSE)
  }
  roles <- tibble::tibble(
    column = c("rds_uri", "summary_parquet_uri", "lbf_parquet_uri", "full_parquet_uri",
               "covariate_explained_variants_uri", "variant_position_summary_uri",
               "outcome_index_uri", "log_uri", "provenance_uri"),
    suffix = c("_susie.rds", ".parquet", ".lbf_variable.parquet", ".full_susie.parquet",
               ".covariate_explained_variants.tsv", ".variant_position_summary.tsv",
               ".outcome_index.tsv", ".log", ".provenance.json"),
    group = c("rds", "parquet", "parquet", "parquet", "tsv", "tsv", "tsv", "log", "provenance"))
  uri_groups <- list(rds = gene_rds_uris, parquet = gene_parquet_uris,
                     tsv = gene_tsv_uris, log = gene_log_uris, provenance = gene_provenance_uris)
  for (group in names(uri_groups)) {
    uris <- uri_groups[[group]]
    if (!is.character(uris) || anyNA(uris) || any(!nzchar(uris)) || any(grepl("[\t\r\n]", uris))) {
      stop("Artifact URI metadata must contain nonempty single-line strings for group: ", group, call. = FALSE)
    }
    # basename() parses text only. URI metadata must never reach filesystem helpers.
    supplied_basenames <- basename(uris)
    if (anyDuplicated(supplied_basenames)) {
      stop("Duplicate artifact URI role in group: ", group, call. = FALSE)
    }
    expected_basenames <- unlist(lapply(roles$suffix[roles$group == group], function(suffix) {
      paste0(genes$gene_key, suffix)
    }), use.names = FALSE)
    if (any(!supplied_basenames %in% expected_basenames)) {
      stop("Foreign artifact or incorrect URI metadata group: ", group, call. = FALSE)
    }
    if (!setequal(supplied_basenames, expected_basenames)) {
      stop("Missing artifact URI roles in group: ", group, call. = FALSE)
    }
  }
  for (i in seq_len(nrow(roles))) {
    uris <- uri_groups[[roles$group[[i]]]]
    positions <- match(paste0(genes$gene_key, roles$suffix[[i]]), basename(uris))
    genes[[roles$column[[i]]]] <- uris[positions]
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  chromosomes <- sort(unique(as.integer(genes$chromosome)))
  files <- stats::setNames(file.path(output_dir, paste0("chr", chromosomes, ".gene_outputs.tsv")),
                          as.character(chromosomes))
  for (chromosome in names(files)) {
    readr::write_tsv(dplyr::filter(genes, .data$chromosome == .env$chromosome),
                     files[[chromosome]], quote = "needed")
  }
  message("Wrote gene-output manifests for ", length(files), " chromosomes and ", nrow(genes), " genes")
  files
}

# Combine native full tables one gene and one Arrow record batch at a time.
# The index supplies ownership; output URI metadata is never opened here.
aggregatePgenFullSusie <- function(outcome_index_file, full_susie_files, output_file,
                                  batch_size = 16384L) {
  if (!is.numeric(batch_size) || length(batch_size) != 1L || is.na(batch_size) ||
      !is.finite(batch_size) || batch_size < 1 || batch_size != floor(batch_size) ||
      batch_size > .Machine$integer.max) {
    stop("batch_size must be one positive integer", call. = FALSE)
  }
  batch_size <- as.integer(batch_size)
  index <- pgenReadTextTable(outcome_index_file,
    numeric_columns = "credible_set_count", logical_columns = "converged")
  owner_columns <- c("batch_id", "gene_id", "outcome_id", "phenotype_id", "gene_key")
  pgenRequireColumns(index, c(owner_columns, "chromosome", "status", "converged",
    "credible_set_count"), "Full-SuSiE outcome index")
  if (!nrow(index)) stop("Full-SuSiE outcome index must not be empty", call. = FALSE)
  for (column in c(owner_columns, "chromosome")) pgenCheckIds(index[[column]], paste("Outcome index", column))
  pgenAutosome(index$chromosome)
  if (anyDuplicated(pgenOwnershipKey(index)) || anyDuplicated(index$outcome_id) ||
      anyDuplicated(index$phenotype_id)) {
    stop("Full-SuSiE outcome index contains duplicate outcome IDs", call. = FALSE)
  }
  if (anyNA(index$status) || any(!index$status %in% c("CONVERGED", "NONCONVERGED")) ||
      anyNA(index$converged) || any(index$converged != (index$status == "CONVERGED"))) {
    stop("Full-SuSiE outcome index has invalid status or converged values", call. = FALSE)
  }
  if (any(!is.finite(index$credible_set_count)) || any(index$credible_set_count < 0) ||
      any(index$credible_set_count != floor(index$credible_set_count))) {
    stop("Full-SuSiE outcome index has invalid credible_set_count", call. = FALSE)
  }
  gene_keys <- vapply(index$gene_id, function(gene) paste0("gene-",
    digest::digest(charToRaw(enc2utf8(gene)), algo = "sha256", serialize = FALSE)),
    character(1), USE.NAMES = FALSE)
  if (any(index$gene_key != gene_keys)) {
    stop("Outcome index gene_key does not match the original gene_id hash", call. = FALSE)
  }
  genes <- index |>
    dplyr::mutate(chromosome = as.character(pgenAutosome(chromosome))) |>
    dplyr::distinct(gene_id, gene_key, batch_id, chromosome)
  if (anyDuplicated(genes$gene_id) || anyDuplicated(genes$gene_key)) {
    stop("Each gene must have exactly one chromosome and batch", call. = FALSE)
  }
  if (!is.character(full_susie_files) || anyNA(full_susie_files) ||
      anyDuplicated(full_susie_files)) {
    stop("Full-SuSiE inputs must be unique local files; duplicate paths are not allowed", call. = FALSE)
  }
  for (path in full_susie_files) pgenLocalFile(path, "Full-SuSiE input")
  normalized_inputs <- normalizePath(full_susie_files, mustWork = TRUE)
  if (anyDuplicated(normalized_inputs) || anyDuplicated(basename(full_susie_files))) {
    stop("Duplicate localized full-SuSiE files or gene basenames", call. = FALSE)
  }
  expected_names <- paste0(genes$gene_key, ".full_susie.parquet")
  if (length(full_susie_files) != nrow(genes) ||
      !setequal(basename(full_susie_files), expected_names)) {
    stop("Supply exactly one full-SuSiE file for each gene; missing or foreign gene basename", call. = FALSE)
  }
  full_susie_files <- full_susie_files[match(expected_names, basename(full_susie_files))]
  if (!is.character(output_file) || length(output_file) != 1L || is.na(output_file) ||
      !nzchar(output_file) || grepl("[\t\r\n]", output_file) ||
      grepl("^[A-Za-z][A-Za-z0-9+.-]*://", output_file)) {
    stop("Full-SuSiE output must be one local path", call. = FALSE)
  }
  output_path <- normalizePath(output_file, mustWork = FALSE)
  if (output_path %in% c(normalized_inputs, normalizePath(outcome_index_file, mustWork = TRUE))) {
    stop("Full-SuSiE output must not replace an input file", call. = FALSE)
  }
  text_columns <- c("molecular_trait_id", "variant", "chromosome", "ref", "alt",
    "cs_id", "cs_index", "region")
  numeric_columns <- c("pip", "z", "posterior_mean", "posterior_sd",
    "X_column_scale_factors", paste0("alpha", seq_len(10L)),
    paste0("mu_", seq_len(10L)), paste0("mu2_", seq_len(10L)))
  required <- c(text_columns, "position", "low_purity", numeric_columns)
  numeric_types <- c("double", "float", "int8", "int16", "int32", "int64",
    "uint8", "uint16", "uint32", "uint64")
  inspect_file <- function(path) {
    reader <- arrow::ParquetFileReader$create(path, mmap = FALSE)
    schema <- reader$GetSchema()
    names <- schema$names
    if (anyDuplicated(names) || !all(required %in% names) || any(owner_columns %in% names)) {
      stop("Full-SuSiE schema has missing, duplicate, or reserved ownership columns", call. = FALSE)
    }
    types <- vapply(schema$fields, function(field) field$type$ToString(), character(1))
    names(types) <- names
    if (any(types[text_columns] != "string") || any(!types[numeric_columns] %in% numeric_types) ||
        !types[["position"]] %in% c("string", numeric_types) ||
        !types[["low_purity"]] %in% if (reader$num_rows) "bool" else c("bool", "string")) {
      stop("Full-SuSiE schema has incompatible column types", call. = FALSE)
    }
    # Normalize only optional additions; all other columns retain exact checks.
    if ("phenotype_residual_SD" %in% names && !types[["phenotype_residual_SD"]] %in% numeric_types) {
      stop("Full-SuSiE requires numeric phenotype_residual_SD", call. = FALSE)
    }
    if ("feature_converged" %in% names && types[["feature_converged"]] != "bool") {
      stop("Full-SuSiE requires logical feature_converged", call. = FALSE)
    }
    optional <- c("phenotype_residual_SD", "feature_converged")
    schema <- arrow::schema(c(schema$fields[!names %in% optional],
      list(arrow::field("phenotype_residual_SD", arrow::float64()),
           arrow::field("feature_converged", arrow::boolean()))))
    types <- c(types[!names(types) %in% optional], phenotype_residual_SD = "double", feature_converged = "bool")
    list(schema = schema, types = types, rows = reader$num_rows)
  }
  # Only footer schemas and row counts are retained across genes.
  descriptions <- lapply(full_susie_files, inspect_file)
  populated <- which(vapply(descriptions, function(item) item$rows > 0, logical(1)))
  prototype <- descriptions[[if (length(populated)) populated[[1L]] else 1L]]
  for (description in descriptions) {
    if (!identical(description$schema$names, prototype$schema$names)) {
      stop("Full-SuSiE files have incompatible schema columns or column order", call. = FALSE)
    }
    different <- names(description$types)[description$types != prototype$types]
    # The legacy zero-row constructor uses integer positions and text purity.
    compatible_empty <- !description$rows &&
      all(different %in% c("position", "low_purity")) &&
      (!"position" %in% different ||
        (description$types[["position"]] == "int32" && prototype$types[["position"]] == "string") ||
        (description$types[["position"]] == "string" && prototype$types[["position"]] == "int32"))
    if (length(different) && !compatible_empty) {
      stop("Full-SuSiE files have incompatible schema types", call. = FALSE)
    }
  }
  output_schema <- arrow::schema(c(prototype$schema$fields,
    lapply(owner_columns, function(column) arrow::field(column, arrow::utf8()))))
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  temporary_output <- tempfile(".full_susie-", tmpdir = dirname(output_file), fileext = ".parquet.tmp")
  sink <- writer <- batch_reader <- NULL
  complete <- FALSE
  on.exit({
    if (!is.null(batch_reader)) try(batch_reader$Close(), silent = TRUE)
    if (!is.null(writer)) try(writer$Close(), silent = TRUE)
    if (!is.null(sink)) try(sink$close(), silent = TRUE)
    if (!complete && file.exists(temporary_output)) unlink(temporary_output)
  }, add = TRUE)
  sink <- arrow::FileOutputStream$create(temporary_output)
  writer <- arrow::ParquetFileWriter$create(output_schema, sink,
    properties = arrow::ParquetWriterProperties$create(output_schema$names))
  total_rows <- 0
  for (i in seq_along(full_susie_files)) {
    own_index <- dplyr::filter(index, gene_key == genes$gene_key[[i]])
    seen_outcomes <- character()
    # Cross-batch duplicate keys are limited to the current gene, then released.
    seen_variants <- new.env(hash = TRUE, parent = emptyenv())
    seen_residual_sd <- new.env(hash = TRUE, parent = emptyenv())
    dataset <- arrow::open_dataset(full_susie_files[[i]], format = "parquet")
    batch_reader <- arrow::Scanner$create(dataset, use_threads = FALSE,
      batch_size = batch_size,
      fragment_scan_options = arrow::ParquetFragmentScanOptions$create(pre_buffer = FALSE))$ToRecordBatchReader()
    repeat {
      batch <- batch_reader$read_next_batch()
      if (is.null(batch)) break
      if (!batch$num_rows) next
      table <- tibble::as_tibble(as.data.frame(batch))
      pgenCheckIds(table$molecular_trait_id, "Full-SuSiE molecular_trait_id")
      owner <- match(table$molecular_trait_id, own_index$outcome_id)
      if (anyNA(owner)) stop("Full-SuSiE rows contain foreign gene outcome ownership", call. = FALSE)
      table <- pgenCheckPhenotypeResidualSD(table, "Full-SuSiE", seen_residual_sd)
      table <- pgenCheckFeatureConverged(table, own_index$converged[owner], "Full-SuSiE")
      for (column in setdiff(text_columns, c("cs_id", "cs_index"))) {
        pgenCheckIds(table[[column]], paste("Full-SuSiE", column))
      }
      for (column in c("cs_id", "cs_index")) {
        pgenCheckIds(table[[column]][!is.na(table[[column]])], paste("Full-SuSiE", column))
      }
      if (any(!grepl("^[^_]+_[0-9]+_[^_]+_.+$", table$variant))) {
        stop("Full-SuSiE has malformed variant IDs", call. = FALSE)
      }
      # The native formatter uses the third and fourth underscore-delimited
      # fields. Additional ID suffix fields do not change the exported allele.
      variant_fields <- stringr::str_split_fixed(table$variant, stringr::fixed("_"), 5L)
      if (any(table$ref != variant_fields[, 3L])) {
        stop("Full-SuSiE ref allele does not match the variant ID", call. = FALSE)
      }
      if (any(table$alt != variant_fields[, 4L])) {
        stop("Full-SuSiE alt allele does not match the variant ID", call. = FALSE)
      }
      variant_chromosome <- sub("^(?i:chr)+", "", sub("_.*$", "", table$variant), perl = TRUE)
      table_chromosome <- sub("^(?i:chr)+", "", table$chromosome, perl = TRUE)
      if (any(pgenAutosome(table_chromosome) != pgenAutosome(variant_chromosome))) {
        stop("Full-SuSiE chromosome does not match the variant ID", call. = FALSE)
      }
      position <- suppressWarnings(as.numeric(table$position))
      variant_position <- suppressWarnings(as.numeric(sub("^[^_]+_([0-9]+)_.*$", "\\1", table$variant)))
      if (any(!is.finite(position)) || any(position < 1) || any(position != floor(position)) ||
          any(position != variant_position)) {
        stop("Full-SuSiE position does not match the variant ID", call. = FALSE)
      }
      if (!all(vapply(table[numeric_columns], is.numeric, logical(1))) || !is.logical(table$low_purity)) {
        stop("Full-SuSiE rows require numeric fields and logical low_purity", call. = FALSE)
      }
      if (any(!is.finite(table$X_column_scale_factors)) || any(table$X_column_scale_factors <= 0)) {
        stop("Full-SuSiE X_column_scale_factors must be positive finite values", call. = FALSE)
      }
      if (any(!is.finite(table$pip)) || any(table$pip < 0 | table$pip > 1)) {
        stop("Full-SuSiE pip must contain finite values from zero to one", call. = FALSE)
      }
      row_keys <- paste0(nchar(table$molecular_trait_id, type = "bytes"), ":",
        table$molecular_trait_id, nchar(table$variant, type = "bytes"), ":", table$variant)
      if (anyDuplicated(row_keys) || any(vapply(row_keys, exists, logical(1),
          envir = seen_variants, inherits = FALSE))) {
        stop("Full-SuSiE contains duplicate outcome variant rows", call. = FALSE)
      }
      for (row_key in row_keys) assign(row_key, TRUE, envir = seen_variants)
      seen_outcomes <- union(seen_outcomes, table$molecular_trait_id)
      for (column in owner_columns) table[[column]] <- own_index[[column]][owner]
      table <- table[output_schema$names]
      writer$WriteBatch(arrow::RecordBatch$create(table, schema = output_schema), chunk_size = batch_size)
      total_rows <- total_rows + nrow(table)
    }
    batch_reader$Close()
    batch_reader <- NULL
    if (any(!own_index$outcome_id[own_index$credible_set_count > 0] %in% seen_outcomes)) {
      stop("Full-SuSiE is missing rows for an outcome with credible sets", call. = FALSE)
    }
    rm(seen_variants, seen_residual_sd, dataset)
    message("Aggregated full-SuSiE gene ", genes$gene_id[[i]])
  }
  if (total_rows == 0) {
    writer$WriteTable(arrow::Table$create(schema = output_schema), chunk_size = batch_size)
  }
  writer$Close()
  writer <- NULL
  sink$close()
  sink <- NULL
  if (!file.rename(temporary_output, output_file)) {
    stop("Cannot move the complete full-SuSiE Parquet file to its output path", call. = FALSE)
  }
  complete <- TRUE
  message("Aggregated ", total_rows, " full-SuSiE rows from ", nrow(genes), " genes")
  invisible(output_file)
}
