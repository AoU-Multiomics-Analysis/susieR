run_memory_debug_tests <- function(repo_dir = ".") {
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(tidyr))
  suppressPackageStartupMessages(library(purrr))
  functions_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH", unset = file.path(repo_dir, "R", "utils"))
  core <- new.env(parent = environment())
  for (filename in c("SusiePerformance.R", "ImportFunctions.R", "SusieFunctions.R", "OptParser.R")) {
    source(file.path(functions_path, filename), local = core)
  }
  diagnostics <- new.env(parent = core)
  source(file.path(functions_path, "MemoryHelpers.R"), local = diagnostics)
  original_options <- options("susier.memory_debug")
  on.exit(options(original_options), add = TRUE)

  set.seed(613L)
  n <- 120L
  p <- 20L
  sample_ids <- paste0("sample", seq_len(n))
  genotype <- matrix(rbinom(n * p, 2L, 0.3), nrow = p,
                     dimnames = list(paste0("chr1_", 1000000L + seq_len(p), "_A_G"), sample_ids))
  sex <- rbinom(n, 1L, 0.5)
  covariates <- cbind(sex = sex, pc = rnorm(n))
  rownames(covariates) <- sample_ids
  phenotype <- genotype[8L, ] + rnorm(n, sd = 0.5)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = rbind(gene1 = phenotype)),
    colData = data.frame(sample_id = sample_ids, genotype_id = sample_ids,
                         qtl_group = "ALL", row.names = sample_ids),
    rowData = data.frame(phenotype_id = "gene1", chromosome = "1",
                         phenotype_pos = 1000010L, row.names = "gene1"))
  explained_genotype <- genotype
  explained_genotype[1L, ] <- sex

  # Removing either gate must fail here before a costly reporter can run.
  # Keep the real fit and instrument only reporting and explicit collections.
  run_fit <- function(debug, genotypes, trap_reporters = FALSE) {
    options(susier.memory_debug = debug)
    scope <- new.env(parent = if (trap_reporters) core else diagnostics)
    if (trap_reporters) {
      scope$debug_list_locals <- function(...) stop("Default path called a costly memory reporter")
      scope$debug_print_specific <- scope$debug_list_locals
    }
    gc_calls <- 0L
    scope$gc <- function(...) {
      gc_calls <<- gc_calls + 1L
      base::gc(...)
    }
    fit_function <- core$finemapPhenotype
    environment(fit_function) <- scope
    messages <- character()
    set.seed(521L)
    output <- capture.output(
      fit <- withCallingHandlers(
        fit_function("gene1", se, genotypes, covariates, 1000L),
        message = function(condition) {
          messages <<- c(messages, conditionMessage(condition))
          invokeRestart("muffleMessage")
        }))
    list(fit = fit, rng = .Random.seed, messages = messages,
         output = output, gc_calls = gc_calls)
  }
  progress <- c("Processing phenotype:", "Using pre-loaded genotype matrix",
                "Residualizing gene expression", "Subsetting genotype matrix",
                "Standardizing genotypes", "Computing gt_hat")
  for (genotypes in list(genotype, explained_genotype)) {
    default <- run_fit(NULL, genotypes, trap_reporters = TRUE)
    disabled <- run_fit(FALSE, genotypes)
    enabled <- run_fit(TRUE, genotypes)
    stopifnot("Debugging must not change the complete fitted object" =
                identical(default$fit, enabled$fit) && identical(disabled$fit, enabled$fit),
              "Debugging must not change RNG state" =
                identical(default$rng, enabled$rng) && identical(disabled$rng, enabled$rng),
              "Default execution must retain the three cleanup collections" =
                identical(default$gc_calls, 3L) && identical(disabled$gc_calls, 3L),
              "Debug execution must add the two diagnostic collections" =
                identical(enabled$gc_calls, 5L),
              "Detailed memory messages must be off by default" =
                !any(grepl("DEBUG:|Estimated gt_hat memory:|GC before", default$messages)),
              "Explicitly disabling debugging must suppress detailed messages" =
                !any(grepl("DEBUG:|Estimated gt_hat memory:|GC before", disabled$messages)))
    for (prefix in progress) {
      stopifnot("Ordinary progress messages must remain available" =
                  any(startsWith(default$messages, prefix)))
    }
    for (prefix in c("DEBUG: before imputation", "DEBUG: before computing gt_hat",
                     "Estimated gt_hat memory:", "GC before imputation:", "GC before gt_hat:",
                     "---- gt_matrix ----", "---- gt_std ----")) {
      stopifnot("Opt-in must run and report the detailed memory diagnostics" =
                  any(startsWith(enabled$messages, prefix)))
    }
    stopifnot("Opt-in must print object sizes and collection results" =
                any(grepl("size_mb", enabled$output)), any(grepl("Ncells", enabled$output)))
    expected_explained <- if (identical(genotypes, explained_genotype)) rownames(genotypes)[1L] else character()
    stopifnot("The fixtures must exercise both ordinary and legacy fallback fits" =
                identical(enabled$fit$covariate_explained_variants$variant_id, expected_explained))
  }

  parser <- optparse::OptionParser(option_list = core$Optlist())
  stopifnot("Legacy CLI memory diagnostics must default to off" =
              identical(optparse::parse_args(parser, args = character())$memory_debug, FALSE),
            "Legacy CLI must accept an explicit opt-in" =
              identical(optparse::parse_args(parser, args = c("--memory_debug", "true"))$memory_debug, TRUE),
            "Legacy CLI must accept an explicit opt-out" =
              identical(optparse::parse_args(parser, args = c("--memory_debug", "false"))$memory_debug, FALSE))
  message("Memory debug tests passed.")
}
