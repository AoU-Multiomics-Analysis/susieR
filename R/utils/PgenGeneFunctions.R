# Impute after the existing hardcall filters and fit sample selection. ReadList
# uses stored dosages; this path expects integer-valued hardcall dosages.
imputePgenHardcalls <- function(reader, mat, missing_value = -1) {
  if (!inherits(reader, "PgenReader") || isTRUE(reader$closed)) {
    stop("Native imputation requires an open PgenReader", call. = FALSE)
  }
  missing_counts <- matrixStats::rowCounts(mat, value = NA_real_) +
    matrixStats::rowCounts(mat, value = missing_value, na.rm = TRUE)
  # An entirely missing variant has no mean and must not reach ReadList().
  all_missing <- which(missing_counts == ncol(mat))
  if (length(all_missing)) {
    message("Removing ", length(all_missing), " entirely missing PGEN variants")
    keep <- missing_counts < ncol(mat)
    mat <- mat[keep, , drop = FALSE]
    missing_counts <- missing_counts[keep]
  }
  missing_rows <- which(missing_counts > 0L & missing_counts < ncol(mat))
  if (!length(missing_rows)) return(mat)

  variant_ids <- rownames(mat)[missing_rows]
  metadata_ids <- reader$variants$variant_id
  if (anyDuplicated(metadata_ids[metadata_ids %in% variant_ids])) {
    stop("Native imputation has ambiguous duplicate variant keys", call. = FALSE)
  }
  variant_indices <- match(variant_ids, metadata_ids)
  sample_indices <- match(colnames(mat), reader$raw_sample_ids)
  if (anyNA(variant_indices) || anyNA(sample_indices) ||
      anyDuplicated(colnames(mat)) || !length(sample_indices)) {
    stop("Native imputation sample or variant IDs do not match the PGEN", call. = FALSE)
  }
  sorted_samples <- sort(sample_indices)
  sample_order <- match(sample_indices, sorted_samples)
  pgen <- pgenlibr::NewPgen(reader$pgen_path, pvar = reader$pvar,
    raw_sample_ct = reader$raw_sample_count,
    sample_subset = as.integer(sorted_samples))
  on.exit(pgenlibr::ClosePgen(pgen), add = TRUE)
  message("Mean-imputing ", length(missing_rows), " PGEN variants across ",
          ncol(mat), " selected samples with pgenlibr")
  # Small read blocks limit temporary matrices. No genotypes persist between fits.
  for (start in seq.int(1L, length(missing_rows), by = 256L)) {
    block <- seq.int(start, min(start + 255L, length(missing_rows)))
    native <- pgenlibr::ReadList(pgen, as.integer(variant_indices[block]),
                               meanimpute = TRUE)
    mat[missing_rows[block], ] <- t(native[sample_order, , drop = FALSE])
  }
  mat
}

# Fit a prepared window while keeping phenotype preparation local to the fit.
# Return NULL only when the existing full-window numerical guard is required.
finemapPreparedPgenPhenotype <- function(phenotype_id, se, reader, covariates,
                                         cis_distance, metadata, genotype_cache,
                                         cache_outcome_id, variant_report_file,
                                         timing_callback = NULL) {
  timing_event <- function(stage, event) {
    if (!is.null(timing_callback)) timing_callback(stage, event)
    invisible(NULL)
  }
  timing_stage <- "data_preparation"
  timing_event(timing_stage, "start")
  on.exit({
    if (!is.null(timing_stage)) timing_event(timing_stage, "failed")
  }, add = TRUE)
  cache <- genotype_cache
  if (!is.environment(cache) || !inherits(cache, "PgenGenotypeCache") ||
      isTRUE(cache$cleared) || !identical(cache$reader, reader)) {
    stop("Prepared genotype cache does not match the open PGEN reader", call. = FALSE)
  }
  gene_vector <- eQTLUtils::extractPhentypeFromSE(phenotype_id, se, "counts") %>%
    dplyr::mutate(phenotype_value_std = qnorm((rank(phenotype_value, na.last = "keep") - 0.5) /
                                            sum(!is.na(phenotype_value))))
  sample_ids <- as.character(gene_vector$genotype_id)
  design <- cbind(covariates[sample_ids, , drop = FALSE], 1)
  if (!identical(sample_ids, cache$sample_ids) ||
      !identical(unname(design), unname(cache$covariates_matrix))) {
    stop("Prepared genotype cache sample order or covariates differ from the feature", call. = FALSE)
  }
  position <- as.numeric(metadata$phenotype_pos[[1L]])
  start <- max(1, position - cis_distance)
  end <- position + cis_distance
  expected <- .pgenWindowIndices(reader, metadata$chromosome[[1L]], start, end,
                                 allow_empty = TRUE)
  if (!identical(expected, cache$plan$window_indices[[cache_outcome_id]])) {
    stop("Prepared genotype cache window differs from the feature window", call. = FALSE)
  }
  state <- new.env(parent = emptyenv())
  state$hat <- NULL
  state$legacy_cache <- NULL
  state$covariate_explained_variants <- emptyCovariateExplainedVariants()
  adjustment <- list(covariates = design, projection = cache$projection, state = state)
  expression_vector <- residualizeWithCovariates(adjustment, gene_vector$phenotype_value_std)
  names(expression_vector) <- sample_ids
  if (!is.null(state$hat)) {
    cache$guard_recomputations <- cache$guard_recomputations + 1L
    message("USE existing per-feature numerical guard for phenotype: ", phenotype_id)
    timing_event("data_preparation", "complete")
    timing_stage <- NULL
    return(NULL)
  }
  window <- getPgenGenotypeCacheWindow(cache, cache_outcome_id)
  if (window$guard_required) {
    cache$guard_recomputations <- cache$guard_recomputations + 1L
    message("USE existing per-feature numerical guard for genotypes: ", phenotype_id)
    timing_event("data_preparation", "complete")
    timing_stage <- NULL
    return(NULL)
  }
  if (!length(window$variant_names)) stop("No variant with alt alleles", call. = FALSE)
  message("Using ", cache$mode, " covariate-adjusted genotype cache for ", phenotype_id,
          ": ", length(window$variant_names), " variants")
  reportVariantDirectionCounts(window$variant_names, metadata,
    context = "fine-mapping genotype matrix", output_file = variant_report_file,
    requested_window_start = position - cis_distance, requested_window_end = end)
  timing_event("data_preparation", "complete")
  timing_stage <- "fine_mapping"
  timing_event(timing_stage, "start")
  fitted <- fitSusieModel(window$gt_hat, expression_vector, legacy = FALSE,
    covariate_explained_variants = window$covariate_explained_variants)
  timing_event(timing_stage, "complete")
  timing_stage <- NULL
  fitted$variant_id <- window$variant_names
  fitted$cohort_maf <- window$cohort_maf
  fitted
}

# Per-outcome storage adapter. Statistical preparation stays in finemapPhenotype().
finemapPgenPhenotype <- function(phenotype_id, se, reader, covariates,
                                 cis_distance, missing_code,
                                 ancestry_data = NULL, MAF = 0,
                                 variant_list = NULL, additional_genotypes = NULL,
                                 variant_report_file = NULL, timing_callback = NULL,
                                 genotype_cache = NULL, cache_outcome_id = phenotype_id) {
  if (!is.null(timing_callback) && !is.function(timing_callback)) {
    stop("timing_callback must be a function or NULL", call. = FALSE)
  }
  loading <- TRUE
  if (!is.null(timing_callback)) timing_callback("data_loading", "start")
  on.exit({
    if (loading && !is.null(timing_callback)) timing_callback("data_loading", "failed")
  }, add = TRUE)
  .pgenGeneId(phenotype_id, "phenotype_id")
  metadata <- as.data.frame(SummarizedExperiment::rowData(se)) |>
    dplyr::filter(.data$phenotype_id == .env$phenotype_id)
  if (nrow(metadata) != 1L) {
    stop("Expected one metadata row for phenotype: ", phenotype_id, call. = FALSE)
  }
  if (!is.numeric(cis_distance) || length(cis_distance) != 1L ||
      is.na(cis_distance) || !is.finite(cis_distance) || cis_distance < 0 ||
      cis_distance != floor(cis_distance)) {
    stop("cis_distance must be a nonnegative integer", call. = FALSE)
  }
  needs_full_cohort_maf <- is.null(variant_list) && !is.null(ancestry_data) && MAF > 0
  if (needs_full_cohort_maf &&
      (length(reader$sample_ids) != reader$raw_sample_count ||
       !setequal(reader$sample_ids, reader$raw_sample_ids))) {
    stop("Ancestry MAF requires the complete raw PGEN sample cohort", call. = FALSE)
  }
  if (!is.null(genotype_cache)) {
    if (needs_full_cohort_maf || !is.null(variant_list) || !is.null(additional_genotypes)) {
      stop("Prepared caching is unavailable for this genotype preparation path", call. = FALSE)
    }
    if (!is.null(timing_callback)) timing_callback("data_loading", "complete")
    loading <- FALSE
    cached <- finemapPreparedPgenPhenotype(phenotype_id, se, reader, covariates,
      cis_distance, metadata, genotype_cache, cache_outcome_id, variant_report_file,
      timing_callback = timing_callback)
    if (!is.null(cached)) return(cached)
    loading <- TRUE
    if (!is.null(timing_callback)) timing_callback("data_loading", "start")
  }
  position <- as.numeric(metadata$phenotype_pos[[1L]])
  start <- position - cis_distance
  end <- position + cis_distance
  message("Reading PGEN window for ", phenotype_id, ": ",
          metadata$chromosome[[1L]], ":", max(1, start), "-", end)
  native_first <- is.null(additional_genotypes) && !needs_full_cohort_maf
  feature_samples <- if (native_first) {
    as.character(eQTLUtils::extractPhentypeFromSE(phenotype_id, se, "counts")$genotype_id)
  } else NULL
  if (native_first) {
    message("Mean-imputing with pgenlibr on the first read across ",
            length(feature_samples), " selected samples")
  }
  raw <- readPgenWindow(reader, metadata$chromosome[[1L]], max(1, start), end,
                        as_double = TRUE,
                        missing_code = if (is.null(additional_genotypes)) "NA" else missing_code,
                        meanimpute = native_first, sample_ids = feature_samples)
  message("PGEN window dimensions: ", nrow(raw), " variants x ", ncol(raw), " samples")
  filtered <- filterMAF(raw, ancestry_data, MAF_threshold = MAF, variant_list = variant_list)
  message("PGEN variants after MAF/list filtering: ", nrow(filtered))
  if (!is.null(timing_callback)) timing_callback("data_loading", "complete")
  loading <- FALSE
  genotype_imputer <- if (!native_first && is.null(additional_genotypes)) {
    function(mat, missing_value) imputePgenHardcalls(reader, mat, missing_value)
  } else NULL
  finemapPhenotype(phenotype_id, se, filtered, covariates, cis_distance,
                   ancestry_data = ancestry_data, MAF = MAF, variant_list = variant_list,
                   additional_genotypes = additional_genotypes,
                   variant_report_file = variant_report_file, timing_callback = timing_callback,
                   genotype_imputer = genotype_imputer,
                   genotypes_imputed = native_first)
}

.pgenGeneId <- function(value, label) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop(label, " must be one nonempty string", call. = FALSE)
  }
  enc2utf8(value)
}

geneArtifactKey <- function(gene_id) {
  gene_id <- .pgenGeneId(gene_id, "gene_id")
  paste0("gene-", digest::digest(charToRaw(gene_id), algo = "sha256", serialize = FALSE))
}

# v1 encodes each UTF-8 tuple member as <byte count>:<value>. The first seven
# SHA256 hex digits give a nonnegative integer below 2^28, accepted by set.seed().
# Batch identity, outcome order and output paths are deliberately absent.
geneOutcomeSeed <- function(base_seed, gene_id, outcome_id) {
  if (!is.numeric(base_seed) || length(base_seed) != 1L || is.na(base_seed) ||
      !is.finite(base_seed) || base_seed < 0 || base_seed > .Machine$integer.max ||
      base_seed != floor(base_seed)) {
    stop("base_seed must be one nonnegative R integer", call. = FALSE)
  }
  values <- c("susier-gene-outcome-seed-v1", as.character(as.integer(base_seed)),
              .pgenGeneId(gene_id, "gene_id"), .pgenGeneId(outcome_id, "outcome_id"))
  tuple <- paste0(nchar(values, type = "bytes"), ":", values, collapse = "")
  hash <- digest::digest(charToRaw(tuple), algo = "sha256", serialize = FALSE)
  strtoi(substr(hash, 1L, 7L), base = 16L)
}

.pgenGeneWithSeed <- function(seed, code) {
  previous_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(previous_kind))
    if (had_seed) {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  # Fix all RNG algorithms so the seed has the same meaning for every caller.
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  force(code)
}

# Use one connection so native stdout and R conditions keep their order.
# Warnings remain original conditions; the caller still controls their handling.
.pgenGeneCaptureLog <- function(log_file, code, warning_callback = NULL) {
  output_depth <- sink.number(type = "output")
  connection <- file(log_file, open = "at")
  on.exit({
    while (sink.number(type = "output") > output_depth) sink(type = "output")
    flush(connection)
    close(connection)
  }, add = TRUE)
  sink(connection, type = "output")
  withCallingHandlers(force(code),
    message = function(condition) {
      text <- conditionMessage(condition)
      cat(text, if (endsWith(text, "\n")) "" else "\n", file = connection, sep = "")
      tryInvokeRestart("muffleMessage")
    },
    warning = function(condition) {
      cat("WARNING: ", conditionMessage(condition), "\n", file = connection, sep = "")
      if (!is.null(warning_callback)) warning_callback(condition)
    })
}

.pgenGeneEmptyIndex <- function() {
  tibble::tibble(gene_id = character(), gene_key = character(), outcome_id = character(),
                phenotype_id = character(), seed = integer(), status = character(),
                converged = logical(), credible_set_count = integer(), reason = character(),
                data_loading_seconds = double(), data_preparation_seconds = double(),
                fine_mapping_seconds = double(), outcome_total_seconds = double())
}

.pgenGeneEmptyPositions <- function() {
  tibble::tibble(
    phenotype_id = character(), context = character(), feature_chromosome = character(),
    feature_position = double(), requested_window_start = double(),
    requested_window_end = double(), requested_window_size_bp = double(),
    total_variants_tested = integer(), same_chromosome_parseable_variants = integer(),
    upstream_variants = integer(), downstream_variants = integer(),
    variants_at_feature_coordinate = integer(), uncounted_variants = integer(),
    min_abs_distance_bp = double(), max_abs_distance_bp = double(),
    min_signed_distance_bp = double(), max_signed_distance_bp = double())
}

# Save large fitted objects once, then return only compact batch-level tables.
runPgenGene <- function(gene_id, outcomes, reader, output_root, base_seed,
                        missing_code, provenance = list(), genotype_cache = NULL) {
  gene_started <- proc.time()[["elapsed"]]
  gene_key <- geneArtifactKey(gene_id)
  outcome_ids <- names(outcomes)
  if (!is.list(outcomes) || !length(outcomes) || is.null(outcome_ids) ||
      anyNA(outcome_ids) || any(!nzchar(outcome_ids)) || anyDuplicated(outcome_ids)) {
    stop("outcomes must be a nonempty named list with unique outcome IDs", call. = FALSE)
  }
  if (!is.character(missing_code) || length(missing_code) != 1L ||
      is.na(missing_code) || !missing_code %in% c("NA", "-1")) {
    stop("missing_code must be \"NA\" or \"-1\"", call. = FALSE)
  }
  if (!is.list(provenance)) stop("provenance must be a list", call. = FALSE)
  seeds <- vapply(outcome_ids, function(id) geneOutcomeSeed(base_seed, gene_id, id), integer(1))
  for (id in outcome_ids) {
    outcome <- outcomes[[id]]
    if (!is.list(outcome) ||
        !all(c("phenotype_id", "se", "covariates", "cis_distance", "region") %in% names(outcome))) {
      stop("Outcome ", id, " lacks required inputs", call. = FALSE)
    }
    .pgenGeneId(outcome$phenotype_id, "phenotype_id")
    .pgenGeneId(outcome$region, "region")
  }
  directory <- file.path(output_root, "genes", gene_key)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(directory, gene_key)
  extra_files <- c(variant_positions = paste0(prefix, ".variant_position_summary.tsv"),
                   index = paste0(prefix, ".outcome_index.tsv"),
                   provenance = paste0(prefix, ".provenance.json"),
                   log = paste0(prefix, ".log"))
  # A failed retry must not leave a successful manifest from an earlier attempt.
  if (file.exists(extra_files[["provenance"]])) unlink(extra_files[["provenance"]])
  writeLines(character(), extra_files[["log"]])
  log_line <- function(text, progress = FALSE) {
    cat(text, "\n", file = extra_files[["log"]], append = TRUE, sep = "")
    if (progress) message(text)
  }
  index <- .pgenGeneEmptyIndex()
  gene_complete <- FALSE
  timing_logged <- FALSE
  log_gene_timing <- function() {
    timing <- list(
      data_loading_seconds = sum(index$data_loading_seconds),
      data_preparation_seconds = sum(index$data_preparation_seconds),
      fine_mapping_seconds = sum(index$fine_mapping_seconds),
      gene_total_seconds = proc.time()[["elapsed"]] - gene_started,
      scope = "gene processing and output writing; excludes final hashing and provenance writing")
    log_line(paste0("Gene timing: ", gene_id, "; status: ",
      if (gene_complete) "COMPLETE" else "FAILED", "; ",
      paste(names(timing)[1:4], format(unlist(timing[1:4]), digits = 6L),
            sep = "=", collapse = "; "), "; ", timing$scope))
    timing_logged <<- TRUE
    timing
  }
  on.exit({
    if (!timing_logged) log_gene_timing()
  }, add = TRUE)
  readr::write_tsv(index, extra_files[["index"]], quote = "needed")
  readr::write_tsv(.pgenGeneEmptyPositions(), extra_files[["variant_positions"]],
                   quote = "needed")
  empty_tables <- formatSusieOutputs(list(), tibble::tibble(phenotype_id = character(), region = character()))
  readr::write_tsv(empty_tables$covariate_explained_variants,
                   paste0(prefix, ".covariate_explained_variants.tsv"), quote = "needed")
  results <- list()
  extracted_results <- list()
  log_line(paste0("[GENE] ", gene_id, "; fits=", length(outcomes),
                  "; log=", extra_files[["log"]]), progress = TRUE)
  for (outcome_number in seq_along(outcome_ids)) {
    id <- outcome_ids[[outcome_number]]
    outcome <- outcomes[[id]]
    seed <- seeds[[id]]
    outcome_started <- proc.time()[["elapsed"]]
    stages <- c("data_loading", "data_preparation", "fine_mapping")
    stage_started <- stats::setNames(rep(NA_real_, length(stages)), stages)
    stage_seconds <- stats::setNames(rep(NA_real_, length(stages)), paste0(stages, "_seconds"))
    timing_callback <- function(stage, event) {
      now <- proc.time()[["elapsed"]]
      if (event == "start") {
        stage_started[[stage]] <<- now
      } else {
        seconds <- now - stage_started[[stage]]
        field <- paste0(stage, "_seconds")
        previous <- stage_seconds[[field]]
        stage_seconds[[field]] <<- if (is.finite(previous)) previous + seconds else seconds
        message("Outcome timing: ", id, "; stage: ", stage, "; status: ", event,
                "; elapsed_seconds=", format(seconds, digits = 6L))
      }
      invisible(NULL)
    }
    timing_row <- function() tibble::as_tibble_row(c(as.list(stage_seconds),
      list(outcome_total_seconds = proc.time()[["elapsed"]] - outcome_started)))
    warning_count <- 0L
    fit_progress <- function(status, count = NULL, fitted = NULL, reason = NULL) {
      details <- paste0("[FIT ", outcome_number, "/", length(outcomes), "] ", status, " ", id)
      if (status == "START") {
        if (outcome$phenotype_id != id) {
          details <- paste0(details, "; phenotype=", outcome$phenotype_id)
        }
      } else {
        if (!is.null(count)) details <- paste0(details, "; cs=", count)
        if (!is.null(fitted$n) && length(fitted$n) == 1L) {
          details <- paste0(details, "; n=", fitted$n)
        }
        if (!is.null(fitted$variant_id)) {
          details <- paste0(details, "; p=", length(fitted$variant_id))
        }
        if (!is.null(fitted$niter) && length(fitted$niter) == 1L) {
          details <- paste0(details, "; iter=", fitted$niter)
        }
        seconds <- c(stage_seconds, proc.time()[["elapsed"]] - outcome_started)
        details <- paste0(details, "; ",
          paste(c("load", "prep", "fit", "total"), format(seconds, digits = 3L, trim = TRUE),
                sep = "=", collapse = " s; "), " s; warnings=", warning_count)
        if (!is.null(reason)) details <- paste0(details, "; reason: ", reason)
      }
      log_line(details, progress = TRUE)
    }
    fit_progress("START")
    log_line(paste0("Outcome start: ", id, "; phenotype: ", outcome$phenotype_id, "; seed: ", seed))
    fitted <- tryCatch(
      .pgenGeneCaptureLog(extra_files[["log"]],
        .pgenGeneWithSeed(seed, finemapPgenPhenotype(
          outcome$phenotype_id, outcome$se, reader, outcome$covariates,
          outcome$cis_distance, missing_code, ancestry_data = outcome$ancestry_data,
          MAF = if (is.null(outcome$MAF)) 0 else outcome$MAF,
          variant_list = outcome$variant_list, additional_genotypes = outcome$additional_genotypes,
          variant_report_file = extra_files[["variant_positions"]], timing_callback = timing_callback,
          genotype_cache = genotype_cache, cache_outcome_id = id)),
        warning_callback = function(condition) warning_count <<- warning_count + 1L),
      error = function(condition) condition)
    if (inherits(fitted, "error")) {
      reason <- conditionMessage(fitted)
      index <- dplyr::bind_rows(index, dplyr::bind_cols(tibble::tibble(
        gene_id = gene_id, gene_key = gene_key, outcome_id = id,
        phenotype_id = outcome$phenotype_id, seed = seed, status = "FAILED",
        converged = NA, credible_set_count = NA_integer_, reason = reason), timing_row()))
      readr::write_tsv(index, extra_files[["index"]], quote = "needed")
      fit_progress("FAILED", reason = reason)
      stop("Gene ", gene_id, ", outcome ", id, " failed: ", reason, call. = FALSE)
    }
    results[[id]] <- fitted
    converged <- isTRUE(fitted$converged)
    status <- if (converged) "CONVERGED" else "NONCONVERGED"
    extracted_results[[id]] <- extractResults(fitted)
    exported_sets <- extracted_results[[id]]$cs_df
    count <- if (is.null(exported_sets)) 0L else as.integer(nrow(exported_sets))
    index <- dplyr::bind_rows(index, dplyr::bind_cols(tibble::tibble(
      gene_id = gene_id, gene_key = gene_key, outcome_id = id,
      phenotype_id = outcome$phenotype_id, seed = seed, status = status,
      converged = converged, credible_set_count = count,
      reason = if (converged) "" else "SuSiE did not converge"), timing_row()))
    readr::write_tsv(index, extra_files[["index"]], quote = "needed")
    fit_progress(status, count = count, fitted = fitted)
  }
  regions <- tibble::tibble(
    phenotype_id = outcome_ids,
    region = vapply(outcomes, function(outcome) outcome$region, character(1)))
  tables <- formatSusieOutputs(results, regions, extracted_results = extracted_results)
  files <- c(writeSusieOutputs(results, tables, prefix), extra_files)
  log_line(paste0("Gene complete: ", gene_id))
  gene_complete <- TRUE
  gene_timing <- log_gene_timing()
  hashed_files <- files[names(files) != "provenance"]
  artifacts <- tibble::tibble(
    name = names(hashed_files), basename = basename(hashed_files),
    sha256 = vapply(hashed_files, digest::digest, character(1), algo = "sha256", file = TRUE))
  record <- utils::modifyList(provenance, list(
    schema_version = 1L, gene_id = gene_id, gene_key = gene_key,
    base_seed = as.integer(base_seed), missing_code = missing_code,
    seed_algorithm = "susier-gene-outcome-seed-v1:UTF8-byte-length:SHA256-first7hex",
    rng_kind = c("Mersenne-Twister", "Inversion", "Rejection"),
    outcomes = index, timing = gene_timing, artifacts = artifacts))
  jsonlite::write_json(record, files[["provenance"]], auto_unbox = TRUE, pretty = TRUE, na = "null")
  list(summary = tables$in_cs_variant_df, credible_sets = tables$cs_df,
       diagnostics = tables$covariate_explained_variants, index = index, files = files)
}
