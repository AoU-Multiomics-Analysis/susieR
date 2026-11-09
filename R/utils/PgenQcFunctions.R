# Shared hardcall QC for a complete input PGEN cohort.
# PLINK reads the native genotype records. R streams metadata and count rows.

validatePgenQcThresholds <- function(mac, min_call_rate) {
  scalar_number <- function(value) {
    is.numeric(value) && length(value) == 1L && !is.na(value) && is.finite(value)
  }
  if (!scalar_number(mac) || mac < 0 || mac != floor(mac)) {
    stop("mac must be one finite nonnegative integer", call. = FALSE)
  }
  if (!scalar_number(min_call_rate) || min_call_rate < 0 || min_call_rate > 1) {
    stop("min_call_rate must be one finite number from 0 to 1", call. = FALSE)
  }
  invisible(list(mac = as.numeric(mac), min_call_rate = as.numeric(min_call_rate)))
}

.pgenQcLocalFile <- function(path, label) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop(label, " must be one readable local file", call. = FALSE)
  }
  if (grepl("^[A-Za-z][A-Za-z0-9+.-]*://", path)) {
    stop(label, " localization error: unresolved URI: ", path, call. = FALSE)
  }
  if (!file.exists(path) || dir.exists(path) || file.access(path, 4L) != 0L) {
    stop(label, " must be one readable local file: ", path, call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

.pgenQcFileRecord <- function(path) {
  list(present = TRUE, basename = basename(path),
       size_bytes = unname(file.info(path)$size),
       sha256 = digest::digest(file = path, algo = "sha256", serialize = FALSE))
}

.pgenQcPvarHeader <- function(connection) {
  header <- character()
  repeat {
    line <- readLines(connection, n = 1L, warn = FALSE)
    if (!length(line)) stop("PVAR requires a #CHROM header", call. = FALSE)
    header <- c(header, line)
    if (startsWith(line, "##")) next
    columns <- strsplit(sub("^#", "", trimws(line)), "[[:space:]]+")[[1L]]
    if (!startsWith(line, "#") || length(columns) < 5L ||
        !identical(columns[1:5], c("CHROM", "POS", "ID", "REF", "ALT")) ||
        anyDuplicated(columns)) {
      stop("PVAR requires CHROM, POS, ID, REF and ALT as its first five columns", call. = FALSE)
    }
    return(header)
  }
}

.pgenQcPvarFields <- function(lines) {
  fields <- strsplit(trimws(lines), "[[:space:]]+")
  if (any(lengths(fields) < 5L) || any(startsWith(lines, "#"))) {
    stop("PVAR has invalid variant records", call. = FALSE)
  }
  fields
}

.pgenQcRunPlink <- function(executable, args, log_path) {
  status <- system2(executable, vapply(args, shQuote, character(1)),
                    stdout = log_path, stderr = log_path)
  if (status != 0L) {
    details <- paste(utils::tail(readLines(log_path, warn = FALSE), 20L), collapse = "\n")
    stop("PLINK failed; see ", log_path, "\n", details, call. = FALSE)
  }
  invisible(NULL)
}

filterPgenFile <- function(opt) {
  if (is.null(opt$mac)) opt$mac <- 0
  if (is.null(opt$min_call_rate)) opt$min_call_rate <- .95
  if (is.null(opt$plink2)) opt$plink2 <- "plink2"
  if (is.null(opt$plink_memory_mb)) opt$plink_memory_mb <- 1024
  thresholds <- validatePgenQcThresholds(opt$mac, opt$min_call_rate)
  memory <- opt$plink_memory_mb
  if (!is.numeric(memory) || length(memory) != 1L || is.na(memory) ||
      !is.finite(memory) || memory < 640 || memory != floor(memory)) {
    stop("plink_memory_mb must be one integer of at least 640 MiB", call. = FALSE)
  }
  input <- lapply(c("pgen", "pvar", "psam"), function(name) .pgenQcLocalFile(opt[[name]], toupper(name)))
  names(input) <- c("pgen", "pvar", "psam")
  if (grepl("\\.zst$", input$pvar, ignore.case = TRUE)) {
    stop("Zstandard PVAR files are not supported; supply a plain-text PVAR", call. = FALSE)
  }
  if (!is.character(opt$output_dir) || length(opt$output_dir) != 1L ||
      is.na(opt$output_dir) || !nzchar(opt$output_dir) ||
      grepl("^[A-Za-z][A-Za-z0-9+.-]*://", opt$output_dir)) {
    stop("output_dir must be one local directory", call. = FALSE)
  }
  if (!is.character(opt$plink2) || length(opt$plink2) != 1L ||
      is.na(opt$plink2) || !nzchar(opt$plink2)) {
    stop("plink2 must name an executable", call. = FALSE)
  }
  executable <- if (grepl("/", opt$plink2, fixed = TRUE)) opt$plink2 else Sys.which(opt$plink2)
  if (!nzchar(executable) || !file.exists(executable) || file.access(executable, 1L) != 0L) {
    stop("plink2 must name an executable", call. = FALSE)
  }
  executable <- normalizePath(executable, mustWork = TRUE)
  required <- c("readr", "tibble", "digest", "jsonlite")
  if (any(!vapply(required, requireNamespace, logical(1), quietly = TRUE))) {
    stop("Shared PGEN QC requires readr, tibble, digest and jsonlite", call. = FALSE)
  }
  if (!exists(".pgenReadMetadata", mode = "function")) {
    stop("Source PgenFunctions.R before PgenQcFunctions.R", call. = FALSE)
  }
  samples <- .pgenReadMetadata(input$psam, "IID", "PSAM")
  if (!nrow(samples) || anyNA(samples$IID) || any(!nzchar(samples$IID)) || anyDuplicated(samples$IID)) {
    stop("PSAM requires nonempty unique sample IIDs", call. = FALSE)
  }
  sample_count <- nrow(samples)
  if (!dir.exists(opt$output_dir) && !dir.create(opt$output_dir, recursive = TRUE)) {
    stop("Cannot create output_dir", call. = FALSE)
  }
  output_dir <- normalizePath(opt$output_dir, mustWork = TRUE)
  output <- list(pgen = file.path(output_dir, "filtered.pgen"),
                 pvar = file.path(output_dir, "filtered.pvar"),
                 psam = file.path(output_dir, "filtered.psam"),
                 counts = file.path(output_dir, "variant_qc.tsv.gz"),
                 report = file.path(output_dir, "qc.json"))
  if (any(unlist(input, use.names = FALSE) %in% unlist(output, use.names = FALSE))) {
    stop("QC outputs must not replace input files", call. = FALSE)
  }
  work_dir <- tempfile(".pgen-qc-", tmpdir = output_dir)
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
  scratch_pvar <- file.path(work_dir, "row_ids.pvar")
  kept_pvar <- file.path(work_dir, "retained.pvar")
  keep_ids <- file.path(work_dir, "keep_ids.txt")
  row_id <- function(index) paste0("qc_row_", format(index, scientific = FALSE, trim = TRUE))
  chunk_size <- 1000L

  message("Preparing row IDs for shared PGEN hardcall QC")
  variant_count <- local({
    source <- gzfile(input$pvar, "rt")
    destination <- file(scratch_pvar, "wt")
    on.exit(close(source), add = TRUE)
    on.exit(close(destination), add = TRUE)
    writeLines(.pgenQcPvarHeader(source), destination)
    count <- 0
    repeat {
      lines <- readLines(source, n = chunk_size, warn = FALSE)
      if (!length(lines)) break
      .pgenQcPvarFields(lines)
      indices <- count + seq_along(lines)
      modified <- vapply(seq_along(lines), function(index) {
        sub("^([[:space:]]*[^[:space:]]+[[:space:]]+[^[:space:]]+[[:space:]]+)[^[:space:]]+",
            paste0("\\1", row_id(indices[[index]])), lines[[index]])
      }, character(1))
      writeLines(modified, destination)
      count <- count + length(lines)
    }
    count
  })
  if (variant_count == 0) stop("PVAR has no variant records", call. = FALSE)
  version <- system2(executable, "--version", stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(version, "status"))) stop("Cannot read PLINK version", call. = FALSE)
  version <- paste(version, collapse = "\n")
  common <- c("--pgen", input$pgen, "--pvar", scratch_pvar, "--psam", input$psam,
              "--allow-extra-chr", "--threads", "1", "--memory", as.character(memory))
  counts_prefix <- file.path(output_dir, "plinkcounts")
  message("Counting PGEN hardcalls across all input samples, including nonfounders")
  .pgenQcRunPlink(executable,
    c(common, "--geno-counts", "cols=homref,refalt1,homalt1,missing,nobs", "--out", counts_prefix),
    paste0(counts_prefix, ".stdout.log"))

  message("Applying MAC and call-rate thresholds to supported autosomal records")
  retained <- local({
    source <- gzfile(input$pvar, "rt")
    counts <- file(paste0(counts_prefix, ".gcount"), "rt")
    restored <- file(kept_pvar, "wt")
    selection <- file(keep_ids, "wt")
    report <- gzfile(output$counts, "wb")
    on.exit(close(source), add = TRUE)
    on.exit(close(counts), add = TRUE)
    on.exit(close(restored), add = TRUE)
    on.exit(close(selection), add = TRUE)
    on.exit(close(report), add = TRUE)
    writeLines(.pgenQcPvarHeader(source), restored)
    header <- strsplit(readLines(counts, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]]
    expected <- c("#ID", "HOM_REF_CT", "HET_REF_ALT1_CT", "HOM_ALT1_CT", "MISSING_CT", "OBS_CT")
    if (!identical(header, expected)) stop("PLINK genotype count header is invalid", call. = FALSE)
    columns <- c("variant_index", "chromosome", "position", "variant_id", "AC", "AN", "MAC",
                 "CALL_RATE", "PASS_MAC", "PASS_CALL_RATE", "KEEP", "qc_applied", "reason")
    writeLines(paste(columns, collapse = "\t"), report)
    before <- 0
    after <- 0
    supported_after <- 0
    repeat {
      lines <- readLines(source, n = chunk_size, warn = FALSE)
      if (!length(lines)) break
      fields <- .pgenQcPvarFields(lines)
      count_lines <- readLines(counts, n = length(lines), warn = FALSE)
      if (length(count_lines) != length(lines)) stop("PLINK count rows do not match PVAR", call. = FALSE)
      count_fields <- strsplit(count_lines, "\t", fixed = TRUE)
      indices <- before + seq_along(lines)
      if (any(lengths(count_fields) != length(expected)) ||
          !identical(vapply(count_fields, `[[`, character(1), 1L), row_id(indices))) {
        stop("PLINK count row order does not match PVAR", call. = FALSE)
      }
      chromosome <- vapply(fields, `[[`, character(1), 1L)
      position <- vapply(fields, `[[`, character(1), 2L)
      raw_id <- vapply(fields, `[[`, character(1), 3L)
      alt <- vapply(fields, `[[`, character(1), 5L)
      supported <- .pgenNormalizeChromosome(chromosome) %in% as.character(seq_len(22L)) &
        !grepl(",", alt, fixed = TRUE)
      ac <- an <- mac <- call_rate <- rep(NA_real_, length(lines))
      pass_mac <- pass_rate <- rep(NA, length(lines))
      keep <- rep(TRUE, length(lines))
      for (index in which(supported)) {
        value <- suppressWarnings(as.numeric(count_fields[[index]][-1L]))
        if (anyNA(value) || any(!is.finite(value)) || any(value < 0) || any(value != floor(value)) ||
            sum(value[1:3]) != value[[5L]] || value[[5L]] + value[[4L]] != sample_count) {
          stop("Invalid diploid hardcall counts for variant ", indices[[index]], call. = FALSE)
        }
        an[[index]] <- 2 * value[[5L]]
        ac[[index]] <- value[[2L]] + 2 * value[[3L]]
        mac[[index]] <- min(ac[[index]], an[[index]] - ac[[index]])
        call_rate[[index]] <- value[[5L]] / sample_count
        pass_mac[[index]] <- mac[[index]] >= thresholds$mac
        pass_rate[[index]] <- call_rate[[index]] >= thresholds$min_call_rate
        keep[[index]] <- pass_mac[[index]] && pass_rate[[index]]
      }
      reason <- ifelse(supported, "applied", "unsupported_chromosome_or_multiallelic")
      chunk <- tibble::tibble(variant_index = indices, chromosome = chromosome,
        position = position, variant_id = raw_id, AC = ac, AN = an, MAC = mac,
        CALL_RATE = call_rate, PASS_MAC = pass_mac, PASS_CALL_RATE = pass_rate,
        KEEP = keep, qc_applied = supported, reason = reason)
      readr::write_tsv(chunk, report, col_names = FALSE, append = TRUE, na = "NA")
      writeLines(lines[keep], restored)
      writeLines(row_id(indices[keep]), selection)
      before <- before + length(lines)
      after <- after + sum(keep)
      supported_after <- supported_after + sum(keep & supported)
    }
    if (length(readLines(counts, n = 1L, warn = FALSE)) || before != variant_count) {
      stop("PLINK count rows do not match PVAR", call. = FALSE)
    }
    list(total = after, supported = supported_after)
  })
  if ((thresholds$mac > 0 || thresholds$min_call_rate > 0) && retained$supported == 0) {
    stop("No supported variants remain after shared PGEN QC", call. = FALSE)
  }
  filtered_prefix <- file.path(output_dir, "filtered")
  message("Writing filtered native PGEN records")
  .pgenQcRunPlink(executable,
    c(common, "--extract", keep_ids, "--make-pgen", "--out", filtered_prefix),
    file.path(output_dir, "filter.stdout.log"))
  local({
    filtered <- file(output$pvar, "rt")
    selected <- file(keep_ids, "rt")
    on.exit(close(filtered), add = TRUE)
    on.exit(close(selected), add = TRUE)
    .pgenQcPvarHeader(filtered)
    repeat {
      lines <- readLines(filtered, n = chunk_size, warn = FALSE)
      if (!length(lines)) break
      ids <- vapply(.pgenQcPvarFields(lines), `[[`, character(1), 3L)
      if (!identical(ids, readLines(selected, n = length(ids), warn = FALSE))) {
        stop("Filtered PGEN variant row order does not match the selected records", call. = FALSE)
      }
    }
    if (length(readLines(selected, n = 1L, warn = FALSE))) {
      stop("Filtered PGEN has fewer records than the selected PVAR", call. = FALSE)
    }
  })
  filtered_samples <- .pgenReadMetadata(output$psam, "IID", "Filtered PSAM")
  if (!identical(samples$IID, filtered_samples$IID) ||
      ("FID" %in% names(samples) &&
       !identical(samples$FID, if ("FID" %in% names(filtered_samples)) filtered_samples$FID else rep("0", sample_count)))) {
    stop("Filtered PGEN sample order does not match input PSAM", call. = FALSE)
  }
  if (!file.copy(kept_pvar, output$pvar, overwrite = TRUE) ||
      !file.copy(input$psam, output$psam, overwrite = TRUE)) {
    stop("Cannot restore retained PVAR and PSAM metadata", call. = FALSE)
  }
  report <- list(schema_version = 1L, count_source = "pgen_hardcalls",
    sample_scope = "all_input_pgen_samples", include_nonfounders = TRUE,
    sample_count = sample_count, variant_count_before = variant_count,
    variant_count_after = retained$total, thresholds = thresholds,
    source_files = lapply(input, .pgenQcFileRecord),
    filtered_files = lapply(output[c("pgen", "pvar", "psam")], .pgenQcFileRecord),
    counts_file = .pgenQcFileRecord(output$counts), plink_version = version)
  jsonlite::write_json(report, output$report, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  unlink(paste0(counts_prefix, ".gcount"))
  message("Shared PGEN QC complete: ", retained$total, " of ", variant_count, " records retained")
  invisible(output)
}
