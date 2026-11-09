# CPU limits must bound all threaded tools, including nested cgroup quotas.
benchmark_library <- Sys.getenv("BENCH_LIBRARY", unset = "")
if (nzchar(benchmark_library)) .libPaths(c(benchmark_library, .libPaths()))
args <- commandArgs(trailingOnly = TRUE)
repo_dir <- normalizePath(if (length(args)) args[[1L]] else getwd(), mustWork = TRUE)
source(file.path(repo_dir, "tests", "test_helpers.R"))
runtime_path <- file.path(Sys.getenv("SUSIER_FUNCTIONS_PATH",
  unset = file.path(repo_dir, "R", "utils")), "PgenRuntime.R")
stopifnot("PGEN runtime helper is required" = file.exists(runtime_path))
source(runtime_path)

run_pgen_runtime_tests <- function() {
  directory <- tempfile("pgen runtime '")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  run_named_test("invalid thread requests fail before runtime changes", {
    for (request in list(-1, 1.5, NA_real_, Inf, "2", c(1, 2), NULL)) {
      expect_error_message(configurePgenThreads(request), "threads must be")
    }
  })
  run_named_test("visible cgroup v2 ancestors enforce the lowest CPU quota", {
    mount <- file.path(directory, "unified")
    task <- file.path(mount, "parent", "task")
    dir.create(task, recursive = TRUE)
    writeLines("max 100000", file.path(mount, "cpu.max"))
    writeLines("250000 100000", file.path(mount, "parent", "cpu.max"))
    writeLines("800000 100000", file.path(task, "cpu.max"))
    quotas <- .pgenAncestorCpuQuotas(task, mount, "v2")
    stopifnot(min(quotas) == 2L)
    writeLines("50000 100000", file.path(mount, "parent", "cpu.max"))
    stopifnot(min(.pgenAncestorCpuQuotas(task, mount, "v2")) == 1L)
    writeLines("bad quota", file.path(mount, "parent", "cpu.max"))
    stopifnot(min(.pgenAncestorCpuQuotas(task, mount, "v2")) == 1L)
  })
  run_named_test("cgroup mount roots resolve process paths without reading outside the mount", {
    mount <- file.path(directory, "mounted")
    dir.create(file.path(mount, "task"), recursive = TRUE)
    writeLines("300000 100000", file.path(mount, "cpu.max"))
    writeLines("800000 100000", file.path(mount, "task", "cpu.max"))
    cgroups <- file.path(directory, "cgroup")
    mounts <- file.path(directory, "mountinfo")
    writeLines("0::/container/task", cgroups)
    writeLines(paste("36 25 0:32 /container", gsub(" ", "\\040", mount, fixed = TRUE),
      "rw - cgroup2 cgroup rw"), mounts)
    stopifnot(min(.pgenVisibleCpuQuotas(cgroups, mounts)) == 3L)
    stopifnot(length(.pgenAncestorCpuQuotas(directory, mount, "v2")) == 0L)
    writeLines("0::/container/..", cgroups)
    stopifnot(min(.pgenVisibleCpuQuotas(cgroups, mounts)) == 3L)
  })
  run_named_test("visible cgroup v1 parent quotas also bound the task", {
    mount <- file.path(directory, "cpu")
    task <- file.path(mount, "parent", "task")
    dir.create(task, recursive = TRUE)
    for (path in c(mount, dirname(task), task)) {
      writeLines("100000", file.path(path, "cpu.cfs_period_us"))
    }
    writeLines("-1", file.path(mount, "cpu.cfs_quota_us"))
    writeLines("200000", file.path(dirname(task), "cpu.cfs_quota_us"))
    writeLines("800000", file.path(task, "cpu.cfs_quota_us"))
    cgroups <- file.path(directory, "cgroup-v1")
    mounts <- file.path(directory, "mountinfo-v1")
    writeLines("2:cpu,cpuacct:/parent/task", cgroups)
    writeLines(paste("36 25 0:32 /", gsub(" ", "\\040", mount, fixed = TRUE),
      "rw - cgroup cgroup rw,cpu,cpuacct"), mounts)
    stopifnot(min(.pgenVisibleCpuQuotas(cgroups, mounts)) == 2L)
  })

  # Exercise a fresh R process so initialized native libraries cannot conceal
  # a missing setter. The parent test process keeps its original thread state.
  script <- file.path(directory, "configure.R")
  writeLines(c(
    "args <- commandArgs(trailingOnly = TRUE)",
    "library_path <- Sys.getenv('BENCH_LIBRARY', '')",
    "if (nzchar(library_path)) .libPaths(c(library_path, .libPaths()))",
    "source(args[[1L]])",
    "Sys.unsetenv('OMP_NUM_THREADS')",
    "reference_allowed <- as.integer(parallelly::availableCores())",
    "Sys.setenv(OMP_NUM_THREADS = '1')",
    "settings <- configurePgenThreads(as.numeric(args[[2L]]))",
    "settings$reference_allowed <- reference_allowed",
    "settings$environment <- as.list(Sys.getenv(c('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',",
    "  'MKL_NUM_THREADS', 'VROOM_THREADS', 'OMP_MAX_ACTIVE_LEVELS', 'OMP_NESTED')))",
    "settings$actual_blas <- RhpcBLASctl::blas_get_num_procs()",
    "settings$actual_openmp <- RhpcBLASctl::omp_get_max_threads()",
    "settings$actual_data_table <- data.table::getDTthreads()",
    "settings$actual_arrow <- arrow::cpu_count()",
    "settings$crossproduct <- unname(crossprod(matrix(1:12, nrow = 4L)))[1L]",
    "x <- matrix(sin(seq_len(4000L)), nrow = 200L)",
    "settings$matrix_results <- c(as.numeric(crossprod(x)), as.numeric(x %*% cos(1:20)))",
    "jsonlite::write_json(settings, args[[3L]], auto_unbox = TRUE, na = 'null', digits = NA)"
  ), script)
  invoke <- function(request) {
    result_path <- tempfile(tmpdir = directory, fileext = ".json")
    output <- system2(file.path(R.home("bin"), "Rscript"),
      shQuote(c(script, runtime_path, as.character(request), result_path)),
      env = c("OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1", "MKL_NUM_THREADS=1",
              "VROOM_THREADS=1"), stdout = TRUE, stderr = TRUE)
    stopifnot(is.null(attr(output, "status")))
    value <- jsonlite::read_json(result_path, simplifyVector = TRUE)
    stopifnot(any(grepl("CPU", output, fixed = TRUE)))
    value
  }
  run_named_test("automatic and capped threads configure active libraries in a subprocess", {
    for (request in c(0, 1, .Machine$integer.max)) {
      value <- invoke(request)
      selected <- value$selected_threads
      stopifnot(value$requested_threads == request,
        value$available_threads >= 1L, selected >= 1L,
        selected <= value$available_threads,
        selected == if (request == 0) value$available_threads else
          min(request, value$available_threads),
        value$arrow_threads == selected, value$actual_arrow == selected,
        value$data_table_threads == value$actual_data_table,
        value$data_table_threads <= selected,
        value$environment$OMP_NUM_THREADS == as.character(selected),
        value$environment$OPENBLAS_NUM_THREADS == as.character(selected),
        value$environment$MKL_NUM_THREADS == as.character(selected),
        value$environment$VROOM_THREADS == as.character(selected),
        value$environment$OMP_MAX_ACTIVE_LEVELS == "1",
        value$environment$OMP_NESTED == "FALSE",
        value$crossproduct == 30,
        length(value$cpu_constraints) > 0L,
        is.character(value$blas_library), nzchar(value$blas_library))
      if (!is.null(value$blas_threads)) {
        stopifnot(value$blas_threads == value$actual_blas, value$blas_threads <= selected)
      }
      if (!is.null(value$openmp_threads)) {
        stopifnot(value$openmp_threads == value$actual_openmp,
                  value$openmp_threads <= selected)
      }
    }
  })
  run_named_test("real CPU detection stays within visible task limits", {
    value <- invoke(0)
    constraints <- unlist(value$cpu_constraints, use.names = FALSE)
    stopifnot(value$available_threads <= min(constraints))
    quotas <- .pgenVisibleCpuQuotas()
    if (length(quotas)) stopifnot(value$available_threads <= min(quotas))
    expected <- min(c(value$reference_allowed, quotas))
    stopifnot(value$available_threads == expected)
  })
  run_named_test("single and automatic thread matrix operations are numerically equivalent", {
    single <- invoke(1)
    automatic <- invoke(0)
    stopifnot(isTRUE(all.equal(single$matrix_results, automatic$matrix_results,
                              tolerance = 1e-12)))
  })
  message("All PGEN runtime tests passed.")
}
run_pgen_runtime_tests()
