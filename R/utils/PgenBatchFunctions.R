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

# Build one column prototype without letting zero-row legacy types set the
# types of populated data. bind_rows supplies the common types and missing cells.
.pgenAggregatePrototype <- function(tables, populated = which(vapply(tables, nrow, integer(1)) > 0L)) {
  columns <- unique(unlist(lapply(tables, names), use.names = FALSE))
  prototype <- if (length(populated)) {
    dplyr::bind_rows(lapply(tables[populated], function(table) table[FALSE, ]))
  } else tables[[1L]][FALSE, ]
  for (table in tables) {
    for (column in setdiff(names(table), names(prototype))) prototype[[column]] <- table[[column]][FALSE]
  }
  prototype[columns]
}

.pgenBindAggregateTables <- function(tables) {
  prototype <- .pgenAggregatePrototype(tables)
  tables <- lapply(tables, function(table) {
    if (!nrow(table)) {
      for (column in intersect(names(table), names(prototype))) table[[column]] <- prototype[[column]]
    }
    table
  })
  dplyr::bind_rows(tables)
}

# Index metadata only fills missing cells. Conflicting values have no lookup
# entry; they never select a result row or prevent the files from being merged.
.pgenAggregateLookup <- function(index) {
  columns <- intersect(c("batch_id", "gene_id", "outcome_id", "phenotype_id", "gene_key", "converged"), names(index))
  if (!"outcome_id" %in% names(index)) return(list())
  stats::setNames(lapply(columns, function(column) {
    values <- dplyr::distinct(index[, unique(c("outcome_id", column)), drop = FALSE])
    keep <- !is.na(values$outcome_id) & !duplicated(values$outcome_id) &
      !duplicated(values$outcome_id, fromLast = TRUE)
    stats::setNames(values[[column]][keep], values$outcome_id[keep])
  }), columns)
}

.pgenAddAggregateMetadata <- function(table, lookup, full = FALSE, residual_sd = TRUE) {
  trait <- if ("molecular_trait_id" %in% names(table)) table$molecular_trait_id else {
    if ("outcome_id" %in% names(table)) table$outcome_id else NULL
  }
  if (is.null(trait)) trait <- rep(NA_character_, nrow(table))
  if (residual_sd && !"phenotype_residual_SD" %in% names(table)) {
    table$phenotype_residual_SD <- rep(NA_real_, nrow(table))
  }
  fallback <- unname(lookup$converged[as.character(trait)])
  if (is.null(fallback)) fallback <- rep(NA, nrow(table))
  if (!"feature_converged" %in% names(table)) table$feature_converged <- rep(NA, nrow(table))
  missing <- is.na(table$feature_converged)
  table$feature_converged[missing] <- fallback[missing]
  if (full) {
    for (column in c("batch_id", "gene_id", "outcome_id", "phenotype_id", "gene_key")) {
      if (!column %in% names(table)) {
        value <- unname(lookup[[column]][as.character(trait)])
        table[[column]] <- if (is.null(value)) rep(NA_character_, nrow(table)) else value
      }
    }
  }
  table
}

aggregatePgenBatches <- function(expected_outcomes, summary_files, credible_set_files,
                                index_files, diagnostic_files, output_dir) {
  message("Loading supplied batch tables")
  # Retained as a readable CLI/WDL input for compatibility; it does not filter
  # or validate the supplied result rows.
  pgenLocalFile(expected_outcomes, "Expected outcomes input")
  read_parquet <- function(path) tibble::as_tibble(arrow::read_parquet(pgenLocalFile(path)))
  index <- .pgenBindAggregateTables(lapply(index_files, pgenReadTextTable,
    numeric_columns = c("seed", "credible_set_count", "data_loading_seconds",
      "data_preparation_seconds", "fine_mapping_seconds", "outcome_total_seconds"),
    logical_columns = "converged"))
  lookup <- .pgenAggregateLookup(index)
  summary <- .pgenAddAggregateMetadata(
    .pgenBindAggregateTables(lapply(summary_files, read_parquet)), lookup)
  cs <- .pgenAddAggregateMetadata(
    .pgenBindAggregateTables(lapply(credible_set_files, read_parquet)), lookup, residual_sd = FALSE)
  diagnostic <- .pgenBindAggregateTables(lapply(diagnostic_files, pgenReadTextTable,
    numeric_columns = c("input_scale", "residual_sd", "residual_sd_ratio")))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  files <- c(summary = file.path(output_dir, "susie.parquet"),
             credible_sets = file.path(output_dir, "credible_sets.parquet"),
             index = file.path(output_dir, "outcome_index.tsv"),
             diagnostics = file.path(output_dir, "covariate_explained_variants.tsv"))
  arrow::write_parquet(summary, files[["summary"]])
  arrow::write_parquet(cs, files[["credible_sets"]])
  readr::write_tsv(index, files[["index"]], quote = "needed")
  readr::write_tsv(diagnostic, files[["diagnostics"]], quote = "needed")
  message("Combined ", nrow(summary), " summary rows, ", nrow(cs),
    " credible-set rows, ", nrow(index), " index rows, and ", nrow(diagnostic), " diagnostic rows")
  invisible(files)
}

# Artifact locations are output metadata. Only the combined outcome index is read.
writeChromosomeGeneManifests <- function(outcome_index_file, gene_rds_uris,
                                        gene_parquet_uris, gene_tsv_uris,
                                        gene_log_uris, gene_provenance_uris,
                                        output_dir) {
  index <- pgenReadTextTable(outcome_index_file, logical_columns = "converged")
  for (column in c("chromosome", "gene_id", "gene_key", "batch_id")) {
    if (!column %in% names(index)) index[[column]] <- rep(NA_character_, nrow(index))
  }
  genes <- index %>%
    dplyr::mutate(chromosome = sub("^chr", "", as.character(chromosome)),
      chromosome = dplyr::coalesce(chromosome, "unknown")) %>%
    dplyr::distinct(chromosome, gene_id, gene_key, batch_id)
  if (!nrow(genes)) return(stats::setNames(character(), character()))
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
  for (i in seq_len(nrow(roles))) {
    uris <- as.character(uri_groups[[roles$group[[i]]]])
    locations <- dplyr::distinct(tibble::tibble(filename = basename(uris), uri = uris))
    known <- !duplicated(locations$filename) & !duplicated(locations$filename, fromLast = TRUE)
    locations <- locations[known, ]
    positions <- match(paste0(genes$gene_key, roles$suffix[[i]]), locations$filename)
    genes[[roles$column[[i]]]] <- locations$uri[positions]
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  chromosomes <- sort(unique(genes$chromosome))
  # URI contents remain text. Chromosome labels are encoded only for filenames.
  labels <- vapply(chromosomes, utils::URLencode, character(1), reserved = TRUE, repeated = TRUE)
  files <- stats::setNames(file.path(output_dir, paste0("chr", labels, ".gene_outputs.tsv")), chromosomes)
  for (i in seq_along(files)) {
    chromosome <- chromosomes[[i]]
    readr::write_tsv(dplyr::filter(genes, .data$chromosome == .env$chromosome),
                     files[[i]], quote = "needed")
  }
  message("Wrote gene-output manifests for ", length(files), " chromosomes and ", nrow(genes), " genes")
  files
}

# Combine the supplied full tables in file order, one Arrow batch at a time.
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
  lookup <- .pgenAggregateLookup(index)
  if (!length(full_susie_files)) stop("Supply at least one full-SuSiE input file", call. = FALSE)
  for (path in full_susie_files) pgenLocalFile(path, "Full-SuSiE input")
  if (!is.character(output_file) || length(output_file) != 1L || is.na(output_file) ||
      !nzchar(output_file) || grepl("[\t\r\n]", output_file) ||
      grepl("^[A-Za-z][A-Za-z0-9+.-]*://", output_file)) {
    stop("Full-SuSiE output must be one local path", call. = FALSE)
  }
  output_path <- normalizePath(output_file, mustWork = FALSE)
  if (output_path %in% c(normalizePath(full_susie_files, mustWork = TRUE),
      normalizePath(outcome_index_file, mustWork = TRUE))) {
    stop("Full-SuSiE output must not replace an input file", call. = FALSE)
  }
  descriptions <- lapply(full_susie_files, function(path) {
    reader <- arrow::ParquetFileReader$create(path, mmap = FALSE)
    empty <- tibble::as_tibble(as.data.frame(arrow::Table$create(schema = reader$GetSchema())))
    list(rows = reader$num_rows, table = empty)
  })
  populated <- which(vapply(descriptions, function(item) item$rows > 0, logical(1)))
  source_prototype <- .pgenAggregatePrototype(lapply(descriptions, `[[`, "table"), populated)
  prototype <- .pgenAddAggregateMetadata(source_prototype, lookup, full = TRUE)
  output_schema <- arrow::Table$create(prototype)$schema
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
  for (path in full_susie_files) {
    dataset <- arrow::open_dataset(path, format = "parquet")
    batch_reader <- arrow::Scanner$create(dataset, use_threads = FALSE,
      batch_size = batch_size,
      fragment_scan_options = arrow::ParquetFragmentScanOptions$create(pre_buffer = FALSE))$ToRecordBatchReader()
    repeat {
      batch <- batch_reader$read_next_batch()
      if (is.null(batch)) break
      if (!batch$num_rows) next
      table <- tibble::as_tibble(as.data.frame(batch))
      for (column in intersect(c("phenotype_residual_SD", "feature_converged"), names(source_prototype))) {
        if (!column %in% names(table)) table[[column]] <- source_prototype[[column]][rep(NA_integer_, nrow(table))]
      }
      table <- .pgenAddAggregateMetadata(table, lookup, full = TRUE)
      table <- dplyr::bind_rows(prototype, table)[output_schema$names]
      writer$WriteBatch(arrow::RecordBatch$create(table, schema = output_schema), chunk_size = batch_size)
      total_rows <- total_rows + nrow(table)
    }
    batch_reader$Close()
    batch_reader <- NULL
    rm(dataset)
    message("Combined full-SuSiE input: ", basename(path))
  }
  if (total_rows == 0) writer$WriteTable(arrow::Table$create(prototype), chunk_size = batch_size)
  writer$Close()
  writer <- NULL
  sink$close()
  sink <- NULL
  if (!file.rename(temporary_output, output_file)) {
    stop("Cannot move the complete full-SuSiE Parquet file to its output path", call. = FALSE)
  }
  complete <- TRUE
  message("Combined ", total_rows, " full-SuSiE rows from ", length(full_susie_files), " input files")
  invisible(output_file)
}
