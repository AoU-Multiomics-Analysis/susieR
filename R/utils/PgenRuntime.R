# Use the CPUs available to this task, including limits on visible parent
# cgroups. Native libraries need runtime setters as well as child-process env.
.pgenRuntimeReadLines <- function(path) {
  tryCatch(suppressWarnings(readLines(path, warn = FALSE)),
           error = function(error) character())
}

.pgenAncestorCpuQuotas <- function(directory, root, version = c("v2", "v1")) {
  version <- match.arg(version)
  root <- normalizePath(root, mustWork = FALSE)
  directory <- normalizePath(directory, mustWork = FALSE)
  prefix <- if (root == "/") "/" else paste0(root, "/")
  if (!identical(directory, root) && !startsWith(directory, prefix)) return(integer())
  limits <- integer()
  repeat {
    path <- file.path(directory, if (version == "v2") "cpu.max" else "cpu.cfs_quota_us")
    if (file.exists(path)) {
      fields <- strsplit(trimws(paste(.pgenRuntimeReadLines(path), collapse = " ")),
                         "[[:space:]]+")[[1L]]
      unlimited <- length(fields) > 0L &&
        fields[[1L]] == if (version == "v2") "max" else "-1"
      if (!unlimited) {
        quota <- suppressWarnings(as.numeric(fields))
        if (version == "v1") {
          quota <- c(quota, suppressWarnings(as.numeric(.pgenRuntimeReadLines(
            file.path(directory, "cpu.cfs_period_us")))))
        }
        # Unreadable or malformed visible limits must not permit host-wide use.
        count <- if (length(quota) == 2L && all(is.finite(quota)) && all(quota > 0)) {
          max(1, floor(quota[[1L]] / quota[[2L]]))
        } else 1L
        limits[[paste0("ancestor:", path)]] <- as.integer(count)
      }
    }
    if (identical(directory, root)) break
    parent <- dirname(directory)
    if (identical(parent, directory)) break
    directory <- parent
  }
  limits
}

.pgenVisibleCpuQuotas <- function(cgroup_file = "/proc/self/cgroup",
                                  mountinfo_file = "/proc/self/mountinfo") {
  groups <- strsplit(.pgenRuntimeReadLines(cgroup_file), ":", fixed = TRUE)
  mounts <- strsplit(.pgenRuntimeReadLines(mountinfo_file), " - ", fixed = TRUE)
  limits <- integer()
  # mountinfo encodes spaces and other special path characters with octal codes.
  decode_path <- function(path) {
    for (pair in list(c("\\040", " "), c("\\011", "\t"),
                      c("\\012", "\n"), c("\\134", "\\"))) {
      path <- gsub(pair[[1L]], pair[[2L]], path, fixed = TRUE)
    }
    path
  }
  for (mount in mounts) {
    if (length(mount) != 2L) next
    before <- strsplit(mount[[1L]], " ", fixed = TRUE)[[1L]]
    after <- strsplit(mount[[2L]], " ", fixed = TRUE)[[1L]]
    if (length(before) < 5L || length(after) < 3L) next
    version <- if (after[[1L]] == "cgroup2") "v2" else
      if (after[[1L]] == "cgroup" &&
          "cpu" %in% strsplit(after[[3L]], ",", fixed = TRUE)[[1L]]) "v1" else NULL
    if (is.null(version)) next
    root <- decode_path(before[[4L]])
    point <- decode_path(before[[5L]])
    for (group in groups) {
      if (length(group) != 3L) next
      matches <- if (version == "v2") identical(group[[2L]], "") else
        "cpu" %in% strsplit(group[[2L]], ",", fixed = TRUE)[[1L]]
      if (!matches) next
      limits <- c(limits, .pgenAncestorCpuQuotas(point, point, version))
      path <- group[[3L]]
      relative <- if (root == "/") sub("^/", "", path) else
        if (path == root || path == "/") "" else
          if (startsWith(path, paste0(root, "/"))) substring(path, nchar(root) + 2L) else ""
      directory <- file.path(point, relative)
      # A cgroup namespace can hide the process path. Its visible mount root
      # still supplies a useful upper bound, even when the child is not visible.
      if (!dir.exists(directory)) directory <- point
      limits <- c(limits, .pgenAncestorCpuQuotas(directory, point, version))
    }
  }
  limits[!duplicated(names(limits))]
}

configurePgenThreads <- function(threads = 0L) {
  if (!is.numeric(threads) || length(threads) != 1L || is.na(threads) ||
      !is.finite(threads) || threads < 0 || threads != floor(threads) ||
      threads > .Machine$integer.max) {
    stop("threads must be one nonnegative integer; 0 selects available CPUs", call. = FALSE)
  }
  for (package in c("parallelly", "RhpcBLASctl", "data.table", "arrow")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("CPU configuration requires package: ", package, call. = FALSE)
    }
  }
  # nproc can honor OMP_NUM_THREADS. The image's initial value of 1 must not
  # become a task CPU limit. Kernel affinity and quota limits remain in force.
  previous_omp <- Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)
  Sys.unsetenv("OMP_NUM_THREADS")
  constraints <- tryCatch(parallelly::availableCores(which = "all"),
    finally = {
      if (is.na(previous_omp)) Sys.unsetenv("OMP_NUM_THREADS") else
        Sys.setenv(OMP_NUM_THREADS = previous_omp)
    })
  constraints <- c(constraints, .pgenVisibleCpuQuotas())
  constraints <- constraints[is.finite(constraints) & constraints >= 1L]
  if (!length(constraints)) constraints <- c(fallback = 1L)
  available <- as.integer(min(constraints))
  selected <- if (threads == 0) available else as.integer(min(threads, available))
  thread_env <- stats::setNames(rep(as.character(selected), 6L),
    c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
      "VECLIB_MAXIMUM_THREADS", "VROOM_THREADS", "BLIS_NUM_THREADS"))
  do.call(Sys.setenv, as.list(thread_env))
  Sys.setenv(OMP_MAX_ACTIVE_LEVELS = "1", OMP_NESTED = "FALSE")
  RhpcBLASctl::blas_set_num_threads(selected)
  RhpcBLASctl::omp_set_num_threads(selected)
  data.table::setDTthreads(selected)
  arrow::set_cpu_count(selected)
  native_count <- function(value) {
    if (length(value) != 1L || is.na(value) || !is.finite(value) || value < 1) NA_integer_ else
      as.integer(value)
  }
  blas_library <- unname(extSoftVersion()["BLAS"])
  if (is.na(blas_library) || !nzchar(blas_library)) blas_library <- "unknown"
  result <- list(requested_threads = as.integer(threads), available_threads = available,
    selected_threads = selected,
    blas_threads = native_count(RhpcBLASctl::blas_get_num_procs()),
    openmp_threads = native_count(RhpcBLASctl::omp_get_max_threads()),
    data_table_threads = as.integer(data.table::getDTthreads()),
    arrow_threads = as.integer(arrow::cpu_count()),
    cpu_constraints = as.list(constraints), blas_library = blas_library)
  message("[RESOURCES] CPU threads: requested=", threads,
          if (threads == 0) " (auto)" else "", "; available=", available,
          "; selected=", selected)
  message("[RESOURCES] Active threads: BLAS=", result$blas_threads,
          "; OpenMP=", result$openmp_threads, "; data.table=", result$data_table_threads,
          "; Arrow=", result$arrow_threads, "; genes run sequentially")
  result
}
