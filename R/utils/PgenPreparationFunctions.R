# Prepare one native PGEN shard for each planned batch from a shared cohort.
readPgenPreparationWorklist <- function(path, cis_distance) {
  pgenLocalFile(path, "Worklist")
  worklist <- readr::read_tsv(path,
    col_types = readr::cols(.default = readr::col_character()), na = character(),
    name_repair = "minimal", show_col_types = FALSE)
  required <- c("batch_id", "chromosome", "gene_id", "outcome_id", "phenotype_id",
    "phenotype_pos", "requested_window_start", "requested_window_end",
    "window_start", "window_end", "region")
  if (!all(required %in% names(worklist)) || anyDuplicated(names(worklist)) ||
      nrow(readr::problems(worklist)) || !nrow(worklist) ||
      anyNA(worklist[required]) || any(!nzchar(as.matrix(worklist[required]))) ||
      anyDuplicated(worklist$outcome_id) || anyDuplicated(worklist$phenotype_id) ||
      dplyr::n_distinct(worklist$batch_id) != 1L ||
      dplyr::n_distinct(worklist$chromosome) != 1L) {
    stop("Invalid batch worklist or duplicate outcome IDs", call. = FALSE)
  }
  pgenAutosome(worklist$chromosome)
  numeric <- c("phenotype_pos", "requested_window_start", "requested_window_end", "window_start", "window_end")
  worklist <- worklist |>
    dplyr::mutate(dplyr::across(dplyr::all_of(numeric), ~ suppressWarnings(as.numeric(.x))))
  values <- as.matrix(worklist[numeric])
  if (any(!is.finite(values)) || any(values != floor(values)) ||
      any(worklist$requested_window_start != worklist$phenotype_pos - cis_distance) ||
      any(worklist$requested_window_end != worklist$phenotype_pos + cis_distance) ||
      any(worklist$window_start != pmax(1, worklist$requested_window_start)) ||
      any(worklist$window_end != worklist$requested_window_end) ||
      any(worklist$window_end < worklist$window_start)) {
    stop("Invalid batch cis windows for the supplied cis distance", call. = FALSE)
  }
  worklist
}

.pgenPreparationWindows <- function(worklist) {
  columns <- c("chromosome", "gene_id", "outcome_id", "phenotype_id", "phenotype_pos",
    "requested_window_start", "requested_window_end", "window_start", "window_end", "region")
  lapply(seq_len(nrow(worklist)), function(index) as.list(worklist[index, columns]))
}

# Read source PVAR once. Open each destination only for its header or a selected chunk.
.pgenPreparationWritePvars <- function(source_path, batches, output_paths, chunk_size = 1000L) {
  source <- gzfile(source_path, "rt")
  on.exit(close(source), add = TRUE)
  header <- .pgenQcPvarHeader(source)
  for (path in output_paths) writeLines(header, path)
  chromosome <- vapply(batches, function(batch) as.character(pgenAutosome(batch$chromosome[[1L]])), character(1))
  lower <- vapply(batches, function(batch) min(batch$window_start), numeric(1))
  upper <- vapply(batches, function(batch) max(batch$window_end), numeric(1))
  counts <- numeric(length(batches))
  repeat {
    lines <- readLines(source, n = chunk_size, warn = FALSE)
    if (!length(lines)) break
    fields <- .pgenQcPvarFields(lines)
    row_chr <- .pgenNormalizeChromosome(vapply(fields, `[[`, character(1), 1L))
    row_pos <- suppressWarnings(as.numeric(vapply(fields, `[[`, character(1), 2L)))
    # Unselected source records can have zero or negative position sentinels.
    if (any(!is.finite(row_pos)) || any(row_pos != floor(row_pos))) {
      stop("Source PVAR has invalid variant positions", call. = FALSE)
    }
    for (chr in intersect(unique(row_chr), chromosome)) {
      rows <- which(row_chr == chr)
      candidates <- which(chromosome == chr & lower <= max(row_pos[rows]) & upper >= min(row_pos[rows]))
      for (index in candidates) {
        batch <- batches[[index]]
        selected <- rep(FALSE, length(rows))
        for (window in seq_len(nrow(batch))) {
          selected <- selected | (row_pos[rows] >= batch$window_start[[window]] &
            row_pos[rows] <= batch$window_end[[window]])
        }
        if (!any(selected)) next
        destination <- file(output_paths[[index]], "at")
        tryCatch(writeLines(lines[rows[selected]], destination), finally = close(destination))
        counts[[index]] <- counts[[index]] + sum(selected)
      }
    }
  }
  counts
}

# Check PLINK row order before restoring the original retained metadata.
.pgenPreparationValidatePvar <- function(actual_path, expected_path, chunk_size = 1000L) {
  actual <- gzfile(actual_path, "rt")
  expected <- gzfile(expected_path, "rt")
  on.exit(close(actual), add = TRUE)
  on.exit(close(expected), add = TRUE)
  .pgenQcPvarHeader(actual)
  .pgenQcPvarHeader(expected)
  repeat {
    observed <- readLines(actual, n = chunk_size, warn = FALSE)
    intended <- readLines(expected, n = chunk_size, warn = FALSE)
    if (length(observed) != length(intended)) {
      stop("Prepared PVAR variant count differs from selected source records", call. = FALSE)
    }
    if (!length(intended)) break
    observed <- .pgenQcPvarFields(observed)
    intended <- .pgenQcPvarFields(intended)
    field <- function(rows, column) vapply(rows, `[[`, character(1), column)
    same <- identical(.pgenNormalizeChromosome(field(observed, 1L)),
      .pgenNormalizeChromosome(field(intended, 1L))) &&
      identical(as.numeric(field(observed, 2L)), as.numeric(field(intended, 2L))) &&
      all(vapply(3:5, function(column) identical(field(observed, column), field(intended, column)), logical(1)))
    if (!same) stop("Prepared PVAR variant row order differs from selected source records", call. = FALSE)
  }
  invisible(TRUE)
}

.pgenPreparationCopyText <- function(source_path, output_path) {
  source <- gzfile(source_path, "rb")
  destination <- file(output_path, "wb")
  on.exit(close(source), add = TRUE)
  on.exit(close(destination), add = TRUE)
  repeat {
    bytes <- readBin(source, "raw", n = 1024L * 1024L)
    if (!length(bytes)) break
    writeBin(bytes, destination)
  }
  invisible(output_path)
}

preparePgenBatches <- function(opt) {
  if (is.null(opt$mac)) opt$mac <- 0
  if (is.null(opt$min_call_rate)) opt$min_call_rate <- .95
  if (is.null(opt$plink2)) opt$plink2 <- "plink2"
  if (is.null(opt$plink_memory_mb)) opt$plink_memory_mb <- 1024
  if (is.null(opt$batch_memory_gb)) opt$batch_memory_gb <- 0L
  if (is.null(opt$genotype_cache_mode)) opt$genotype_cache_mode <- "auto"
  if (is.null(opt$genotype_cache_memory_gb)) opt$genotype_cache_memory_gb <- 0
  if (is.null(opt$cache_ineligible)) opt$cache_ineligible <- FALSE
  if (is.null(opt$has_additional_genotypes)) opt$has_additional_genotypes <- FALSE
  if (!is.logical(opt$has_additional_genotypes) || length(opt$has_additional_genotypes) != 1L ||
      is.na(opt$has_additional_genotypes)) stop("Invalid additional-genotype setting", call. = FALSE)
  if (opt$has_additional_genotypes && identical(as.numeric(opt$batch_memory_gb), 0)) {
    stop("Automatic memory sizing requires a fixed --batch-memory-gb with additional genotypes",
      call. = FALSE)
  }
  if (!is.logical(opt$cache_ineligible) || length(opt$cache_ineligible) != 1L ||
      is.na(opt$cache_ineligible)) stop("Invalid cache eligibility", call. = FALSE)
  validatePgenBatchMemoryOptions(opt$genotype_cache_mode,
    opt$genotype_cache_memory_gb, opt$batch_memory_gb, !opt$cache_ineligible)
  validatePgenQcThresholds(opt$mac, opt$min_call_rate)
  integer <- function(value, lower) is.numeric(value) && length(value) == 1L &&
    is.finite(value) && value >= lower && value <= .Machine$integer.max && value == floor(value)
  if (!integer(opt$cis_distance, 0) || !integer(opt$plink_memory_mb, 640)) {
    stop("Invalid cis distance or PLINK memory (minimum 640 MiB)", call. = FALSE)
  }
  inputs <- stats::setNames(lapply(c("pgen", "pvar", "psam"), function(role)
    .pgenQcLocalFile(opt[[role]], role)), c("pgen", "pvar", "psam"))
  lists <- .pgenQcLocalFile(opt$worklists_file, "Worklists file")
  worklists <- readLines(lists, warn = FALSE)
  if (!length(worklists) || any(!nzchar(worklists)) || anyDuplicated(worklists)) {
    stop("Worklists file must contain distinct nonempty local paths", call. = FALSE)
  }
  worklists <- vapply(worklists, .pgenQcLocalFile, character(1), label = "Worklist")
  batches <- lapply(worklists, readPgenPreparationWorklist, cis_distance = opt$cis_distance)
  ids <- vapply(batches, function(batch) batch$batch_id[[1L]], character(1))
  if (anyDuplicated(ids)) stop("Duplicate batch IDs in preparation worklists", call. = FALSE)
  if (any(!grepl("^[A-Za-z0-9_]+$", ids))) {
    stop("Preparation batch IDs must contain only letters, digits or underscores", call. = FALSE)
  }
  executable <- opt$plink2
  if (!is.character(executable) || length(executable) != 1L || is.na(executable) || !nzchar(executable)) {
    stop("PLINK 2 executable is not available", call. = FALSE)
  }
  if (!grepl("/", executable, fixed = TRUE)) executable <- Sys.which(executable)
  if (!nzchar(executable) || !file.exists(executable) || file.access(executable, 1L) != 0L) {
    stop("PLINK 2 executable is not available", call. = FALSE)
  }
  executable <- normalizePath(executable, mustWork = TRUE)
  output <- opt$output_dir
  if (!is.character(output) || length(output) != 1L || is.na(output) || !nzchar(output) ||
      grepl("[\r\n]", output) || grepl("^[A-Za-z][A-Za-z0-9+.-]*://", output)) {
    stop("output_dir must be one local directory", call. = FALSE)
  }
  if (dir.exists(output) && length(list.files(output, all.files = TRUE, no.. = TRUE))) {
    stop("output_dir must be empty; use a new preparation directory", call. = FALSE)
  }
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  output <- normalizePath(output, mustWork = TRUE)
  batch_dir <- file.path(output, "batches")
  dir.create(batch_dir, showWarnings = FALSE)
  shared_qc <- NULL
  qc_report <- NULL
  source <- inputs
  if (opt$mac > 0 || opt$min_call_rate > 0) {
    qc_opt <- opt
    qc_opt$output_dir <- file.path(output, "qc")
    filtered <- filterPgenFile(qc_opt)
    source <- filtered[c("pgen", "pvar", "psam")]
    qc_report <- filtered$report
    shared_qc <- list(file = .pgenQcFileRecord(qc_report),
      report = jsonlite::read_json(qc_report, simplifyVector = FALSE))
    version <- shared_qc$report$plink_version
  } else {
    message("Shared genotype QC disabled; use original cohort")
    version <- system2(executable, "--version", stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(version, "status"))) stop("Cannot read PLINK version", call. = FALSE)
    version <- paste(version, collapse = "\n")
  }
  # QC already hashed its complete filtered files. Reuse those records.
  source_records <- if (is.null(shared_qc)) lapply(source, .pgenQcFileRecord) else
    shared_qc$report$filtered_files
  source_psam <- .pgenReadMetadata(source$psam, "IID", "PSAM")
  source_samples <- source_psam$IID
  metadata_dir <- tempfile(".pgen-metadata-", tmpdir = output)
  dir.create(metadata_dir)
  on.exit(unlink(metadata_dir, recursive = TRUE), add = TRUE)
  metadata_paths <- file.path(metadata_dir, paste0(ids, ".pvar"))
  message("SELECT original PVAR metadata in one source pass")
  expected_counts <- .pgenPreparationWritePvars(source$pvar, batches, metadata_paths)
  paths <- lapply(c("pgen", "pvar", "psam", "reports"), function(role) character(length(batches)))
  names(paths) <- c("pgen", "pvar", "psam", "reports")
  total_records <- total_bytes <- 0
  memory_reports <- vector("list", length(batches))
  message("PREPARE ", length(batches), " batches from one input cohort")
  for (index in seq_along(batches)) {
    worklist <- batches[[index]]
    prefix <- file.path(batch_dir, ids[[index]])
    regions <- paste0(prefix, ".windows.tsv")
    writePgenBatchRegions(worklist, regions)
    message("EXTRACT ", ids[[index]], ": ", nrow(worklist), " windows; retain all samples")
    .pgenQcRunPlink(executable, c("--pgen", source$pgen, "--pvar", source$pvar,
      "--psam", source$psam, "--extract", "range", regions, "--make-pgen",
      "--threads", "1", "--memory", as.character(opt$plink_memory_mb), "--out", prefix),
      paste0(prefix, ".stdout.log"))
    files <- as.list(stats::setNames(paste0(prefix, c(".pgen", ".pvar", ".psam")), c("pgen", "pvar", "psam")))
    emitted_psam <- .pgenReadMetadata(files$psam, "IID", "Prepared PSAM")
    fid <- function(samples) if ("FID" %in% names(samples)) samples$FID else rep("0", nrow(samples))
    if (!identical(emitted_psam$IID, source_samples) || !identical(fid(emitted_psam), fid(source_psam))) {
      stop("Prepared PGEN IID/FID order differs from source cohort", call. = FALSE)
    }
    .pgenPreparationValidatePvar(files$pvar, metadata_paths[[index]])
    if (!file.copy(metadata_paths[[index]], files$pvar, overwrite = TRUE)) {
      stop("Cannot restore selected source PVAR metadata", call. = FALSE)
    }
    .pgenPreparationCopyText(source$psam, files$psam)
    reader <- openPgenReader(files$pgen, files$pvar, files$psam)
    details <- tryCatch({
      if (!identical(reader$raw_sample_ids, source_samples)) {
        stop("Prepared shard sample order differs from source cohort", call. = FALSE)
      }
      if (nrow(reader$variants) != expected_counts[[index]]) {
        stop("Prepared PGEN variant count differs from selected source records", call. = FALSE)
      }
      # Copy IDs to remove readr growable-vector flags from R serialization.
      memory <- planPgenBatchMemory(reader, worklist, cache_mode = opt$genotype_cache_mode,
        cache_memory_gb = opt$genotype_cache_memory_gb, memory_gb = opt$batch_memory_gb,
        cache_eligible = !opt$cache_ineligible)
      memory_reports[[index]] <- c(list(batch_id = ids[[index]]), memory)
      message("MEMORY ", ids[[index]], ": request=", memory$memory_gb,
        " GB; estimated=", memory$estimated_memory_gb, " GB; samples=", memory$sample_count,
        "; largest_window_variants=", memory$max_window_variants,
        "; cache_mode=", memory$cache_mode)
      if (memory$memory_gb < memory$estimated_memory_gb) {
        message("MEMORY fixed override is below the estimate for ", ids[[index]])
      }
      list(sample_count = reader$raw_sample_count,
           sample_ids_sha256 = digest::digest(as.character(reader$raw_sample_ids)[], algo = "sha256", serializeVersion = 2),
           sample_ids_hash_algorithm = "sha256:R-serialize-v2",
           variant_count = nrow(reader$variants), memory = memory)
    }, finally = closePgenReader(reader))
    report <- c(list(schema_version = 1L, batch_id = ids[[index]],
      chromosome = worklist$chromosome[[1L]], cis_distance = opt$cis_distance,
      worklist = .pgenQcFileRecord(worklists[[index]]), windows = .pgenPreparationWindows(worklist),
      source_filtered_files = source_records, shard_files = lapply(files, .pgenQcFileRecord),
      shared_qc = shared_qc, plink_version = version), details)
    report_path <- paste0(prefix, ".preparation.json")
    jsonlite::write_json(report, report_path, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
    for (role in c("pgen", "pvar", "psam")) paths[[role]][[index]] <- files[[role]]
    paths$reports[[index]] <- report_path
    bytes <- sum(vapply(report$shard_files, function(record) record$size_bytes, numeric(1)))
    total_records <- total_records + details$variant_count
    total_bytes <- total_bytes + bytes
    message("COMPLETE ", ids[[index]], ": ", details$variant_count, " records, ", bytes, " bytes")
  }
  message("COMPLETE preparation: ", length(batches), " batches, ", total_records,
    " extracted records, ", total_bytes, " triplet bytes")
  # Cromwell glob outputs use lexicographic batch-ID order, including when a
  # direct CLI caller supplied its worklists in another order.
  memory <- dplyr::bind_rows(memory_reports[order(ids)])
  memory_gb_file <- file.path(output, "batch_memory_gb.json")
  memory_report <- file.path(output, "batch_memory.tsv")
  jsonlite::write_json(as.integer(memory$memory_gb), memory_gb_file,
    auto_unbox = FALSE, pretty = TRUE)
  readr::write_tsv(memory, memory_report)
  invisible(c(paths, list(whole_pgen = source$pgen, whole_pvar = source$pvar,
    whole_psam = source$psam, qc_report = qc_report,
    memory_gb_file = memory_gb_file, memory_report = memory_report)))
}

# Validate small reports and local shard files before phenotype loading or fits.
readPgenPreparedBatchReport <- function(opt, worklist, reader, inputs = NULL) {
  invalid <- function(detail) stop("Invalid prepared batch report: ", detail, call. = FALSE)
  if (!is.null(opt$mac) && opt$mac > 0) {
    stop("Set --mac to 0 with --prepared-batch-manifest; do not repeat the genotype filter", call. = FALSE)
  }
  path <- pgenLocalFile(opt$prepared_batch_manifest, "Prepared batch report")
  report <- tryCatch(jsonlite::read_json(path, simplifyVector = FALSE),
    error = function(error) invalid(conditionMessage(error)))
  fields <- function(value, required) is.list(value) && !is.null(names(value)) &&
    !anyDuplicated(names(value)) && all(required %in% names(value))
  scalar <- function(value) is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
  count <- function(value) is.numeric(value) && length(value) == 1L && is.finite(value) && value >= 0 && value == floor(value)
  record <- function(value) fields(value, c("present", "basename", "size_bytes", "sha256")) &&
    identical(value$present, TRUE) && scalar(value$basename) &&
    identical(basename(value$basename), value$basename) && count(value$size_bytes) &&
    scalar(value$sha256) && grepl("^[a-fA-F0-9]{64}$", value$sha256)
  same_record <- function(actual, expected, label) {
    if (!record(expected) || !record(actual) ||
        actual$size_bytes != expected$size_bytes || tolower(actual$sha256) != tolower(expected$sha256)) {
      invalid(paste(label, "file hash or size differs"))
    }
  }
  required <- c("schema_version", "batch_id", "chromosome", "cis_distance", "worklist", "windows",
    "source_filtered_files", "shard_files", "shared_qc", "plink_version", "sample_count", "sample_ids_sha256", "sample_ids_hash_algorithm", "variant_count")
  if (!fields(report, required) || !count(report$schema_version) || report$schema_version != 1L) invalid("unsupported or incomplete schema")
  if (!identical(report$batch_id, worklist$batch_id[[1L]]) ||
      !identical(report$chromosome, worklist$chromosome[[1L]]) ||
      !count(report$cis_distance) || report$cis_distance != opt$cis_distance) invalid("batch ID, chromosome or cis distance differs")
  same_record(.pgenQcFileRecord(pgenLocalFile(opt$worklist, "Worklist")), report$worklist, "Worklist")
  expected_windows <- .pgenPreparationWindows(worklist)
  # JSON normalizes integer and double storage. Compare canonical JSON objects.
  canonical <- function(value) jsonlite::toJSON(value, auto_unbox = TRUE, digits = NA)
  if (!identical(canonical(report$windows), canonical(expected_windows))) invalid("cis windows differ")
  for (group in c("source_filtered_files", "shard_files")) {
    if (!fields(report[[group]], c("pgen", "pvar", "psam")) ||
        !all(vapply(report[[group]][c("pgen", "pvar", "psam")], record, logical(1)))) invalid(paste(group, "records are invalid"))
  }
  for (role in c("pgen", "pvar", "psam")) {
    local <- pgenLocalFile(opt[[role]], toupper(role))
    actual <- if (is.null(inputs)) .pgenQcFileRecord(local) else inputs[[role]]
    same_record(actual, report$shard_files[[role]], paste("Shard", role))
  }
  if (!count(report$sample_count) || report$sample_count != reader$raw_sample_count ||
      !identical(report$sample_ids_hash_algorithm, "sha256:R-serialize-v2") ||
      !identical(report$sample_ids_sha256, digest::digest(as.character(reader$raw_sample_ids)[],
        algo = "sha256", serializeVersion = 2))) invalid("full-cohort samples or their order differ")
  if (!count(report$variant_count) || report$variant_count != nrow(reader$variants)) invalid("variant count differs")
  inside <- rep(FALSE, nrow(reader$variants))
  for (index in seq_len(nrow(worklist))) {
    inside <- inside | (reader$variants$chromosome == as.character(pgenAutosome(worklist$chromosome[[index]])) &
      reader$variants$position >= worklist$window_start[[index]] &
      reader$variants$position <= worklist$window_end[[index]])
  }
  if (!all(inside)) invalid("shard contains variants outside the batch windows")
  if (!scalar(report$plink_version)) invalid("PLINK version is missing")
  if (!is.null(report$shared_qc)) {
    if (is.null(opt$genotype_qc_report)) {
      invalid("shared QC requires --genotype-qc-report to verify its file hash and contents")
    }
    if (!fields(report$shared_qc, c("file", "report")) || !record(report$shared_qc$file)) invalid("shared QC binding is invalid")
    qc <- readPgenBatchQcReport(NULL, report$source_filtered_files, 0,
      reader$raw_sample_count, report$shared_qc$report)
    if (report$variant_count > qc$variant_count_after) invalid("shard exceeds shared QC variant count")
    if (!is.null(opt$genotype_qc_report)) {
      qc_path <- pgenLocalFile(opt$genotype_qc_report, "Genotype QC report")
      same_record(.pgenQcFileRecord(qc_path), report$shared_qc$file, "Shared QC report")
      external <- jsonlite::read_json(qc_path, simplifyVector = FALSE)
      if (!identical(canonical(external), canonical(qc))) invalid("embedded shared QC differs from the supplied report")
    }
  } else if (!is.null(opt$genotype_qc_report)) invalid("shared QC is disabled but a QC report was supplied")
  report
}
