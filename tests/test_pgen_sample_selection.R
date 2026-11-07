# Sample selection checks do not run PLINK or fit models.
run_pgen_sample_selection_tests <- function(repo_dir) {
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(tidyr))
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(purrr))
  source(file.path(repo_dir, "tests", "test_helpers.R"), local = TRUE)
  functions_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH",
                               unset = file.path(repo_dir, "R", "utils"))
  previous_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH", unset = NA_character_)
  Sys.setenv(SUSIER_FUNCTIONS_PATH = functions_path)
  on.exit(if (is.na(previous_path)) Sys.unsetenv("SUSIER_FUNCTIONS_PATH") else
    Sys.setenv(SUSIER_FUNCTIONS_PATH = previous_path), add = TRUE)
  scripts_path <- Sys.getenv("SUSIER_SCRIPTS_PATH",
                             unset = file.path(repo_dir, "R", "scripts"))
  source(file.path(scripts_path, "run_pgen_susie_batch.R"), local = TRUE)
  directory <- tempfile("pgen-sample-selection-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  paths <- setNames(file.path(directory, c("phenotypes.bed", "qtl.tsv", "covariates.tsv",
                                         "samples.tsv", "genotypes.pgen")),
                    c("bed", "qtl", "covariates", "samples", "pgen"))
  writeLines(c("#chr\tstart\tend\tphenotype_id\ttrait_only\ts2\tNA\t001\ts1",
               "1\t100\t101\ttrait1\t99\t2\t3\tNA\t1",
               "1\t200\t201\ttrait2\t99\t5\t6\t7\t4"), paths[["bed"]])
  writeLines(c("phenotype_id\tqval", "trait1\t0.01", "trait2\t0.02"), paths[["qtl"]])
  writeLines(c("ID\tcov_extra\tNA\ts2\t001\ts1", "pc1\t9\t3\t2\t1\t0",
               "pc2\t8\t6\t4\t2\t1"), paths[["covariates"]])
  writeLines(c("sample_id", "s1", "001", "s2"), paths[["samples"]])
  writeLines("not opened by the loader", paths[["pgen"]])
  genotype_ids <- c("genotype_only", "s1", "001", "s2", "NA")
  opt <- list(expression_matrix = paths[["bed"]], phenotype_list = paths[["qtl"]],
    covariates = paths[["covariates"]], genotype_matrix = paths[["pgen"]],
    sample_meta = NULL, cv_meta = NULL, out_prefix = "batch", cisdistance = 100L,
    genotype_sample_ids = genotype_ids[-5L])

  run_named_test("inferred samples use the PSAM-order intersection and keep literal IDs", {
    # A phenotype-order or covariate-order intersection must fail this check.
    metadata <- inferPgenSampleMetadata(genotype_ids,
      c("trait_only", "s2", "NA", "001", "s1"), c("cov_extra", "NA", "s2", "001", "s1"))
    expect_identical_value(metadata, tibble::tibble(
      sample_id = c("s1", "001", "s2", "NA"),
      genotype_id = c("s1", "001", "s2", "NA"), qtl_group = "ALL"))
  })
  run_named_test("inference rejects invalid or duplicate IDs before selecting samples", {
    for (ids in list(character(), c("s1", "s1"), c("s1", NA_character_),
                     c("s1", ""), 1:2)) {
      expect_error_message(inferPgenSampleMetadata(ids, c("s1", "s2"), c("s1", "s2")),
                           "Genotype sample IDs")
      expect_error_message(inferPgenSampleMetadata(c("s1", "s2"), ids, c("s1", "s2")),
                           "Molecular trait sample IDs")
    }
    expect_error_message(inferPgenSampleMetadata("s1", "s2", c("s1", "s2")),
                         "No shared genotype and molecular trait sample IDs")
    expect_error_message(inferPgenSampleMetadata(c("s1", "s2"), c("s1", "s2"), "s1"),
                         "Covariates are missing retained sample IDs")
  })
  run_named_test("the loader retains phenotype NA values and matches an explicit intersection", {
    inferred <- LoadData(opt)
    explicit_opt <- opt
    explicit_opt$sample_meta <- paths[["samples"]]
    explicit_opt$genotype_sample_ids <- NULL
    explicit <- LoadData(explicit_opt)
    expect_identical_value(as.data.frame(inferred$sample_metadata),
                           as.data.frame(explicit$sample_metadata), "selected sample metadata")
    expect_identical_value(unname(as.matrix(inferred$expression_matrix)),
                           unname(as.matrix(explicit$expression_matrix)), "molecular trait values")
    expect_identical_value(inferred$covariates_matrix, explicit$covariates_matrix,
                           "covariate values")
    inferred_se <- makePgenBatchExperiment(inferred)
    explicit_se <- makePgenBatchExperiment(explicit)
    expect_identical_value(colnames(inferred_se), c("s1", "001", "s2"))
    expect_identical_value(SummarizedExperiment::assay(inferred_se),
                           SummarizedExperiment::assay(explicit_se), "selected assay")
    expect_identical_value(as.data.frame(SummarizedExperiment::colData(inferred_se)),
                           as.data.frame(SummarizedExperiment::colData(explicit_se)),
                           "selected column metadata")
    stopifnot(is.na(SummarizedExperiment::assay(inferred_se)["trait1", "001"]),
              "cov_extra" %in% rownames(inferred$covariates_matrix))
  })
  run_named_test("the default loader preserves a literal NA sample header", {
    # fread's default missing-value rules otherwise rename this header to V7.
    literal_opt <- opt
    literal_opt$genotype_sample_ids <- c("NA", genotype_ids[-5L])
    data <- LoadData(literal_opt)
    expect_identical_value(data$sample_metadata$sample_id, c("NA", "s1", "001", "s2"))
    se <- makePgenBatchExperiment(data)
    expect_identical_value(colnames(se), c("NA", "s1", "001", "s2"))
    stopifnot(SummarizedExperiment::assay(se)["trait1", "NA"] == 3)
  })
  run_named_test("PGEN metadata keeps parsed IDs without changing sample values", {
    parsed_bed <- file.path(directory, "parsed-phenotypes.bed")
    writeLines(c("#chr\tstart\tend\tphenotype_id\ttrait_only\ts2\tNA\t001\ts1",
                 '" 1 "\t" 100 "\t" 101 "\t" padded "\t99\t2\t3\tNA\t1',
                 '1\t200\t201\t"a""b"\t99\t5\t6\t7\t4',
                 '1\t300\t301\t001\t99\t8\t9\t10\t7',
                 '1\t400\t401\tNA\t99\t11\t12\t13\t10',
                 '1\t500\t501\ttrait1\t99\t14\t15\t16\t13'), parsed_bed)
    parsed_metadata <- tibble::tibble(
      chromosome = rep("1", 5L), start = as.character(seq(100L, 500L, 100L)),
      end = as.character(seq(101L, 501L, 100L)),
      phenotype_id = c("padded", 'a"b', "001", "NA", "trait1"))
    parsed_opt <- opt
    parsed_opt$expression_matrix <- parsed_bed
    parsed_opt$pgen_phenotype_metadata <- readPgenPhenotypeMetadata(parsed_bed)
    expect_identical_value(parsed_opt$pgen_phenotype_metadata, parsed_metadata,
                           "parsed phenotype metadata")
    data <- LoadData(parsed_opt)
    expect_identical_value(data$expression_matrix$phenotype_id,
                           c("padded", 'a"b', "001", "NA", "trait1"))
    expect_identical_value(data$phenotype_meta$phenotype_id,
                           c("padded", 'a"b', "001", "NA", "trait1"))
    expect_identical_value(data$phenotype_meta$phenotype_pos,
                           c(100, 200, 300, 400, 500))
    sample_values <- data$expression_matrix |>
      dplyr::select(-dplyr::all_of(names(data$expression_matrix)[1:4])) |>
      as.matrix()
    expect_identical_value(unname(sample_values), matrix(c(
      99L, 2L, 3L, NA_integer_, 1L,
      99L, 5L, 6L, 7L, 4L,
      99L, 8L, 9L, 10L, 7L,
      99L, 11L, 12L, 13L, 10L,
      99L, 14L, 15L, 16L, 13L), nrow = 5L, byrow = TRUE), "sample values")
    bad_opt <- parsed_opt
    bad_opt$pgen_phenotype_metadata <- dplyr::slice(parsed_metadata, 1:4)
    expect_error_message(LoadData(bad_opt), "PGEN phenotype metadata row count")
    bad_opt$pgen_phenotype_metadata <- dplyr::rename(parsed_metadata, id = phenotype_id)
    expect_error_message(LoadData(bad_opt), "PGEN phenotype metadata schema")
  })
  run_named_test("provided sample metadata retains its order and extra columns", {
    writeLines(c("sample_id\tlabel", "s2\tsecond", "s1\tfirst"), paths[["samples"]])
    explicit_opt <- opt
    explicit_opt$sample_meta <- paths[["samples"]]
    # Internal IDs must have no effect when an explicit list is present.
    explicit_opt$genotype_sample_ids <- c("duplicate", "duplicate")
    data <- LoadData(explicit_opt)
    expect_identical_value(as.data.frame(data$sample_metadata), data.frame(
      sample_id = c("s2", "s1"), label = c("second", "first"),
      genotype_id = c("s2", "s1"), qtl_group = "ALL"))
    expect_identical_value(colnames(makePgenBatchExperiment(data)), c("s2", "s1"))
    legacy_opt <- opt
    legacy_opt$genotype_sample_ids <- NULL
    expect_error_message(LoadData(legacy_opt), "Sample metadata missing")
  })
  run_named_test("inferred samples require covariate coverage without reducing the cohort", {
    bad_opt <- opt
    bad_opt$genotype_sample_ids <- c(genotype_ids, "trait_only")
    expect_error_message(LoadData(bad_opt), "Covariates are missing retained sample IDs")
    bad_opt$genotype_sample_ids <- c("s1", "s1")
    expect_error_message(LoadData(bad_opt), "Genotype sample IDs")
    bad_bed <- file.path(directory, "duplicate-phenotypes.bed")
    lines <- readLines(paths[["bed"]])
    lines[[1L]] <- sub("trait_only", "s2", lines[[1L]], fixed = TRUE)
    writeLines(lines, bad_bed)
    bad_opt <- opt
    bad_opt$expression_matrix <- bad_bed
    expect_error_message(LoadData(bad_opt), "Molecular trait sample IDs")
  })
  run_named_test("inference preserves sample rows with one or zero retained covariates", {
    one_covariate <- file.path(directory, "one-covariate.tsv")
    writeLines(c("ID\tcov_extra\tNA\ts2\t001\ts1", "pc1\t9\t3\t2\t1\t0"), one_covariate)
    one_opt <- opt
    one_opt$covariates <- one_covariate
    data <- LoadData(one_opt)
    expect_identical_value(dim(data$covariates_matrix), c(5L, 1L))
    expect_identical_value(rownames(data$covariates_matrix), c("cov_extra", "NA", "s2", "001", "s1"))
    expect_identical_value(data$sample_metadata$sample_id, c("s1", "001", "s2"))
    writeLines(c("ID\tcov_extra\tNA\ts2\t001\ts1", "pc1\t1\t1\t1\t1\t1",
                 "pc2\t0\t0\t0\t0\t0"), one_covariate)
    intercept_only <- LoadData(one_opt)
    expect_identical_value(dim(intercept_only$covariates_matrix), c(5L, 0L))
    expect_identical_value(rownames(intercept_only$covariates_matrix),
                           c("cov_extra", "NA", "s2", "001", "s1"))
    expect_identical_value(intercept_only$sample_metadata$sample_id, c("s1", "001", "s2"))
  })
  run_named_test("inference rejects duplicate raw covariate headers before row-name repair", {
    duplicate_covariates <- file.path(directory, "duplicate-covariates.tsv")
    writeLines(c("ID\ts1\ts1\t001\ts2\tNA", "pc1\t0\t9\t1\t2\t3",
                 "pc2\t1\t8\t2\t4\t6"), duplicate_covariates)
    legacy_error <- function(validate = NULL) {
      tryCatch({
        if (is.null(validate)) importQtlmapCovariates(duplicate_covariates) else
          importQtlmapCovariates(duplicate_covariates, validate_sample_ids = validate)
        NULL
      }, error = function(condition) list(message = conditionMessage(condition),
                                          class = class(condition)))
    }
    expect_identical_value(legacy_error(FALSE), legacy_error(),
                           "legacy duplicate-header error")
    expect_error_message(importQtlmapCovariates(duplicate_covariates,
                           validate_sample_ids = TRUE), "Covariate sample IDs")
    bad_opt <- opt
    bad_opt$covariates <- duplicate_covariates
    expect_error_message(LoadData(bad_opt), "Covariate sample IDs")
  })
  run_named_test("the batch CLI permits an omitted sample list and validates a supplied path", {
    parser <- optparse::OptionParser(option_list = batch_options)
    parsed <- optparse::parse_args(parser, args = character(), print_help_and_exit = FALSE)
    stopifnot(is.null(parsed$sample_list))
    # Reach worklist validation without PLINK, genotype reads, or model fits.
    batch_opt <- list(pgen = paths[["pgen"]], pvar = paths[["pgen"]], psam = paths[["pgen"]],
      worklist = paths[["qtl"]], phenotype_bed = paths[["bed"]], qtl_results = paths[["qtl"]],
      covariates = paths[["covariates"]], sample_list = NULL, missing_code = "NA",
      output_dir = directory, cis_distance = 100L, base_seed = 1L,
      plink_memory_mb = 1024L, maf = 0)
    expect_error_message(runPgenSusieBatch(batch_opt), "Invalid batch worklist")
    batch_opt$sample_list <- "gs://bucket/unlocalized-samples.tsv"
    expect_error_message(runPgenSusieBatch(batch_opt), "Input localization error")
    batch_opt$sample_list <- file.path(directory, "missing-samples.tsv")
    expect_error_message(runPgenSusieBatch(batch_opt), "Batch inputs must be readable local files")
  })
  invisible(TRUE)
}
