# Test real output connections and a small native fit. No fit output is mocked.
run_pgen_logging_tests <- function(repo_dir, section = c("all", "capture", "fits")) {
  section <- match.arg(section)
  benchmark_library <- Sys.getenv("BENCH_LIBRARY", unset = "")
  if (nzchar(benchmark_library)) .libPaths(c(benchmark_library, .libPaths()))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(tidyr))
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(purrr))
  source(file.path(repo_dir, "tests", "test_helpers.R"), local = TRUE)
  functions_path <- Sys.getenv("SUSIER_FUNCTIONS_PATH",
                               unset = file.path(repo_dir, "R", "utils"))
  for (filename in c("InitFunctions.R", "ImportFunctions.R", "MemoryHelpers.R",
                     "SusiePerformance.R", "SusieFunctions.R", "PgenFunctions.R",
                     "SusieOutputFunctions.R", "PgenGeneFunctions.R")) {
    source(file.path(functions_path, filename), local = TRUE)
  }
  directory <- tempfile("pgen-logging-tests-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)

  if (section != "fits") {
    run_named_test("fit capture keeps output and conditions in emission order", {
      path <- file.path(directory, "ordered.log")
      writeLines("existing", path)
      condition <- structure(list(message = "ordinary warning", call = NULL),
                             class = c("fit_warning", "warning", "condition"))
      warnings <- list()
      warning_count <- 0L
      messages <- character()
      output <- capture.output({
        value <- withCallingHandlers(
          .pgenGeneCaptureLog(path, {
            cat("output before\n")
            message("fit message")
            warning(condition)
            cat("output after\n")
            17L
          }, warning_callback = function(condition) warning_count <<- warning_count + 1L),
          warning = function(condition) {
            warnings[[length(warnings) + 1L]] <<- condition
            invokeRestart("muffleWarning")
          },
          message = function(condition) messages <<- c(messages, conditionMessage(condition)))
      })
      expect_identical_value(output, character(), "captured main stdout")
      expect_identical_value(messages, character(), "captured main messages")
      expect_identical_value(value, 17L, "original return value")
      expect_identical_value(warnings, list(condition), "original warning observed once")
      expect_identical_value(warning_count, 1L, "warning count")
      expect_identical_value(readLines(path), c("existing", "output before", "fit message",
        "WARNING: ordinary warning", "output after"), "ordered append log")
    })
    run_named_test("fit capture preserves warning-to-error behavior and prior output", {
      path <- file.path(directory, "warning-error.log")
      previous <- options(warn = 2L)
      tryCatch({
        expect_error_message(.pgenGeneCaptureLog(path, {
          cat("before warning error\n")
          warning("warning becomes error", call. = FALSE)
          cat("unreachable\n")
        }), "warning becomes error")
      }, finally = options(previous))
      expect_identical_value(readLines(path),
        c("before warning error", "WARNING: warning becomes error"), "failed warning log")
    })
    for (fail in c(FALSE, TRUE)) {
      run_named_test(paste("fit capture restores nested caller sinks after",
                           if (fail) "error" else "success"), {
        check_sinks <- function() {
          before_output <- sink.number(type = "output")
          before_message <- sink.number(type = "message")
          paths <- file.path(directory, paste0(if (fail) "failure" else "success",
                                             c("-outer", "-inner", "-message", "-fit"), ".log"))
          connections <- lapply(paths[1:3], file, open = "wt")
          on.exit({
            while (sink.number(type = "output") > before_output) sink(type = "output")
            if (before_message == 2L) sink(type = "message") else {
              sink(getConnection(before_message), type = "message")
            }
            lapply(connections, close)
          }, add = TRUE)
          sink(connections[[1L]], type = "output")
          sink(connections[[2L]], type = "output")
          sink(connections[[3L]], type = "message")
          expected_output <- sink.number(type = "output")
          expected_message <- sink.number(type = "message")
          fitted <- tryCatch(.pgenGeneCaptureLog(paths[[4L]], {
            cat("native iteration\n")
            message("native detail")
            if (fail) stop("native fit failed", call. = FALSE)
            23L
          }), error = identity)
          stopifnot(sink.number(type = "output") == expected_output,
                    sink.number(type = "message") == expected_message)
          if (fail) stopifnot(inherits(fitted, "error")) else stopifnot(identical(fitted, 23L))
          cat("inner still writable\n")
          message("message still writable")
          sink(type = "output")
          cat("outer still writable\n")
          invisible(lapply(connections, flush))
          expect_identical_value(readLines(paths[[1L]]), "outer still writable")
          expect_identical_value(readLines(paths[[2L]]), "inner still writable")
          expect_identical_value(readLines(paths[[3L]]), "message still writable")
          expect_identical_value(readLines(paths[[4L]]), c("native iteration", "native detail"))
        }
        check_sinks()
      })
    }
  }
  if (section == "capture") return(invisible(TRUE))

  set.seed(907L)
  samples <- paste0("s", seq_len(80L))
  positions <- seq_len(10L) * 10L
  variant_ids <- paste0("chr1_", positions, "_A_G")
  genotypes <- matrix(as.double(rbinom(10L * 80L, 2L, 0.3)), nrow = 10L,
                       dimnames = list(variant_ids, samples))
  paths <- file.path(directory, c("fixture.pgen", "fixture.pvar", "fixture.psam"))
  writer <- pgenlibr::NewPgenWriter(paths[[1L]], ncol(genotypes), nrow(genotypes),
                                   nonref_flags = FALSE)
  tryCatch({
    for (i in seq_len(nrow(genotypes))) pgenlibr::AppendBiallelic(writer, genotypes[i, ])
  }, finally = pgenlibr::ClosePgenWriter(writer))
  writeLines(c("#CHROM\tPOS\tID\tREF\tALT",
               paste("1", positions, variant_ids, "A", "G", sep = "\t")), paths[[2L]])
  writeLines(c("#IID", samples), paths[[3L]])
  reader <- openPgenReader(paths[[1L]], paths[[2L]], paths[[3L]])
  on.exit(closePgenReader(reader), add = TRUE)
  covariates <- matrix(rnorm(80L), ncol = 1L, dimnames = list(samples, "pc"))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = rbind(trait = 3 * genotypes[2L, ] + rnorm(80L, sd = 0.3))),
    colData = data.frame(sample_id = samples, genotype_id = samples,
                         row.names = samples, qtl_group = "ALL"),
    rowData = data.frame(phenotype_id = "trait", chromosome = "1",
                         phenotype_pos = 50L, row.names = "trait"))
  outcomes <- list(expression = list(phenotype_id = "trait", se = se,
    covariates = covariates, cis_distance = 50L, region = "chr1_0_100"))
  # eQTLUtils emits a once-per-session dependency deprecation warning. Warm that
  # path so the fit warning count does not depend on the test runner's order.
  invisible(suppressWarnings(eQTLUtils::extractPhentypeFromSE("trait", se, "counts")))
  run_gene <- function(gene, inputs = outcomes) {
    progress <- character()
    warnings <- list()
    stdout <- capture.output({
      value <- withCallingHandlers(
        tryCatch(runPgenGene(gene, inputs, reader, directory, 13L, "NA"), error = identity),
        message = function(condition) {
          progress <<- c(progress, sub("\n$", "", conditionMessage(condition)))
          invokeRestart("muffleMessage")
        },
        warning = function(condition) {
          warnings[[length(warnings) + 1L]] <<- condition
          invokeRestart("muffleWarning")
        })
    })
    list(value = value, progress = progress, stdout = stdout, warnings = warnings)
  }
  run_named_test("gene console shows progress while the log keeps native fit details", {
    actual <- run_gene("logging-gene")
    stopifnot(!inherits(actual$value, "error"))
    progress <- actual$progress
    expect_identical_value(actual$stdout, character(), "native stdout outside gene log")
    stopifnot(length(progress) == 3L,
              grepl("logging-gene", progress[[1L]], fixed = TRUE),
              grepl("[GENE]", progress[[1L]], fixed = TRUE),
              grepl("fits=1", progress[[1L]], fixed = TRUE),
              grepl(actual$value$files[["log"]], progress[[1L]], fixed = TRUE),
              grepl("[FIT 1/1] START", progress[[2L]], fixed = TRUE),
              grepl("expression", progress[[2L]], fixed = TRUE),
              grepl("trait", progress[[2L]], fixed = TRUE),
              grepl("[FIT 1/1] CONVERGED", progress[[3L]], fixed = TRUE),
              grepl("cs=", progress[[3L]], fixed = TRUE),
              grepl("total=", progress[[3L]], fixed = TRUE),
              grepl("load=", progress[[3L]], fixed = TRUE),
              grepl("prep=", progress[[3L]], fixed = TRUE),
              grepl("fit=", progress[[3L]], fixed = TRUE),
              grepl("warnings=0", progress[[3L]], fixed = TRUE))
    details <- readLines(actual$value$files[["log"]])
    stopifnot(any(grepl("Reading PGEN window", details, fixed = TRUE)),
              any(grepl("Fitting SuSiE model", details, fixed = TRUE)),
              any(grepl("iter.*ELBO|^[[:space:]]*\\[1,\\]", details)),
              any(grepl("Outcome timing:", details, fixed = TRUE)),
              any(grepl("seed:", details, fixed = TRUE)),
              any(grepl("Gene timing:", details, fixed = TRUE)))
    provenance <- jsonlite::read_json(actual$value$files[["provenance"]], simplifyVector = TRUE)
    published_hash <- provenance$artifacts$sha256[provenance$artifacts$name == "log"]
    expect_identical_value(published_hash,
      digest::digest(file = actual$value$files[["log"]], algo = "sha256"), "final gene log hash")
  })
  run_named_test("nonconverged native fits show status and preserve the original warning", {
    # Limit only this test's native fit to one iteration. The real PGEN adapter,
    # native fitting implementation, extraction, and output writing still run.
    production_fit <- runSusieWithFitOptimizations
    tryCatch({
      runSusieWithFitOptimizations <- function(...) production_fit(..., max_iter = 1L)
      actual <- run_gene("nonconverged-logging-gene")
    }, finally = {
      runSusieWithFitOptimizations <- production_fit
    })
    stopifnot(!inherits(actual$value, "error"), length(actual$progress) == 3L,
              grepl("[FIT 1/1] NONCONVERGED", actual$progress[[3L]], fixed = TRUE),
              grepl("expression", actual$progress[[3L]], fixed = TRUE),
              grepl("warnings=1", actual$progress[[3L]], fixed = TRUE),
              length(actual$warnings) == 1L,
              inherits(actual$warnings[[1L]], "warning"),
              identical(actual$value$index$status, "NONCONVERGED"),
              identical(readRDS(actual$value$files[["rds"]])[[1L]]$converged, FALSE))
    warning_text <- conditionMessage(actual$warnings[[1L]])
    details <- readLines(actual$value$files[["log"]])
    # Some native versions also print this warning to stdout. Keep that output;
    # the outer condition handler above must still observe only one warning.
    stopifnot(grepl(paste0("WARNING: ", warning_text), paste(details, collapse = "\n"), fixed = TRUE))
    expect_identical_value(actual$stdout, character(), "nonconverged stdout outside gene log")
  })
  run_named_test("failed fits retain native details and report their outcome context", {
    invalid <- outcomes
    SummarizedExperiment::assay(invalid[[1L]]$se)[1L, ] <- 1
    actual <- run_gene("failed-logging-gene", invalid)
    stopifnot(inherits(actual$value, "error"), length(actual$progress) == 3L,
              grepl("[FIT 1/1] FAILED", actual$progress[[3L]], fixed = TRUE),
              grepl("expression", actual$progress[[3L]], fixed = TRUE),
              grepl("trait", actual$progress[[2L]], fixed = TRUE),
              grepl("total=", actual$progress[[3L]], fixed = TRUE),
              grepl("warnings=", actual$progress[[3L]], fixed = TRUE))
    path <- file.path(directory, "genes", geneArtifactKey("failed-logging-gene"),
                      paste0(geneArtifactKey("failed-logging-gene"), ".log"))
    details <- readLines(path)
    stopifnot(any(grepl("Reading PGEN window", details, fixed = TRUE)),
              any(grepl("Fitting SuSiE model", details, fixed = TRUE)),
              any(grepl("FAILED", details, fixed = TRUE)),
              any(grepl("Gene timing:", details, fixed = TRUE)))
    expect_identical_value(actual$stdout, character(), "failed fit stdout outside gene log")
  })
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  repo_dir <- normalizePath(if (length(args)) args[[1L]] else getwd(), mustWork = TRUE)
  run_pgen_logging_tests(repo_dir)
}
