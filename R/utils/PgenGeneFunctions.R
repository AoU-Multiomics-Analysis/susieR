# Per-outcome storage adapter. Statistical preparation stays in finemapPhenotype().
finemapPgenPhenotype <- function(phenotype_id, se, reader, covariates,
                                 cis_distance, missing_code,
                                 ancestry_data = NULL, MAF = 0,
                                 variant_list = NULL, additional_genotypes = NULL,
                                 variant_report_file = NULL, timing_callback = NULL) {
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
  if (is.null(variant_list) && !is.null(ancestry_data) && MAF > 0 &&
      (length(reader$sample_ids) != reader$raw_sample_count ||
       !setequal(reader$sample_ids, reader$raw_sample_ids))) {
    stop("Ancestry MAF requires the complete raw PGEN sample cohort", call. = FALSE)
  }
  position <- as.numeric(metadata$phenotype_pos[[1L]])
  start <- position - cis_distance
  end <- position + cis_distance
  message("Reading PGEN window for ", phenotype_id, ": ",
          metadata$chromosome[[1L]], ":", max(1, start), "-", end)
  raw <- readPgenWindow(reader, metadata$chromosome[[1L]], max(1, start), end,
                        as_double = TRUE, missing_code = missing_code)
  message("PGEN window dimensions: ", nrow(raw), " variants x ", ncol(raw), " samples")
  filtered <- filterMAF(raw, ancestry_data, MAF_threshold = MAF, variant_list = variant_list)
  message("PGEN variants after MAF/list filtering: ", nrow(filtered))
  if (!is.null(timing_callback)) timing_callback("data_loading", "complete")
  loading <- FALSE
  finemapPhenotype(phenotype_id, se, filtered, covariates, cis_distance,
                   ancestry_data = ancestry_data, MAF = MAF, variant_list = variant_list,
                   additional_genotypes = additional_genotypes,
                   variant_report_file = variant_report_file, timing_callback = timing_callback)
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
                        missing_code, provenance = list()) {
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
  log_line <- function(text) {
    cat(text, "\n", file = extra_files[["log"]], append = TRUE, sep = "")
    message(text)
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
  log_line(paste0("Gene start: ", gene_id, "; outcomes: ", length(outcomes)))
  for (id in outcome_ids) {
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
        stage_seconds[[paste0(stage, "_seconds")]] <<- seconds
        message("Outcome timing: ", id, "; stage: ", stage, "; status: ", event,
                "; elapsed_seconds=", format(seconds, digits = 6L))
      }
      invisible(NULL)
    }
    timing_row <- function() tibble::as_tibble_row(c(as.list(stage_seconds),
      list(outcome_total_seconds = proc.time()[["elapsed"]] - outcome_started)))
    log_line(paste0("Outcome start: ", id, "; phenotype: ", outcome$phenotype_id, "; seed: ", seed))
    fitted <- tryCatch(
      withCallingHandlers(
        .pgenGeneWithSeed(seed, finemapPgenPhenotype(
          outcome$phenotype_id, outcome$se, reader, outcome$covariates,
          outcome$cis_distance, missing_code, ancestry_data = outcome$ancestry_data,
          MAF = if (is.null(outcome$MAF)) 0 else outcome$MAF,
          variant_list = outcome$variant_list, additional_genotypes = outcome$additional_genotypes,
          variant_report_file = extra_files[["variant_positions"]], timing_callback = timing_callback)),
        message = function(condition) cat(conditionMessage(condition), file = extra_files[["log"]], append = TRUE),
        warning = function(condition) cat("WARNING: ", conditionMessage(condition), "\n",
          file = extra_files[["log"]], append = TRUE, sep = "")),
      error = function(condition) condition)
    if (inherits(fitted, "error")) {
      reason <- conditionMessage(fitted)
      index <- dplyr::bind_rows(index, dplyr::bind_cols(tibble::tibble(
        gene_id = gene_id, gene_key = gene_key, outcome_id = id,
        phenotype_id = outcome$phenotype_id, seed = seed, status = "FAILED",
        converged = NA, credible_set_count = NA_integer_, reason = reason), timing_row()))
      readr::write_tsv(index, extra_files[["index"]], quote = "needed")
      log_line(paste0("FAILED outcome ", id, ": ", reason))
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
    log_line(paste0("Outcome complete: ", id, "; ", status, "; credible sets: ", count))
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
