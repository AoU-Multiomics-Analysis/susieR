# ---------- DEBUG HELPERS (put near top of utils/SusieFunctions.R) ----------
# require lobstr for accurate object sizes if available
debug_list_locals <- function(env = parent.frame(), top_n = 50) {
  objs <- ls(envir = env, all.names = TRUE)
  if (length(objs) == 0) {
    message("No local objects in environment")
    return(invisible(data.frame()))
  }
  use_lobstr <- requireNamespace("lobstr", quietly = TRUE)
  info <- lapply(objs, function(nm) {
    obj <- tryCatch(get(nm, envir = env), error = function(e) NULL)
    if (is.null(obj)) return(NULL)
    size_bytes <- if (use_lobstr) as.numeric(lobstr::obj_size(obj)) else as.numeric(object.size(obj))
    cls <- paste(class(obj)[1], collapse = ",")
    dims <- if (!is.null(dim(obj))) paste(dim(obj), collapse = " x ") else length(obj)
    data.frame(name = nm, class = cls, dims = dims, bytes = size_bytes, stringsAsFactors = FALSE)
  })
  info <- do.call(rbind, info[!sapply(info, is.null)])
  info <- info[order(-info$bytes), ]
  info$size_mb <- info$bytes / 1024^2
  print(head(info, n = top_n))
  invisible(info)
}

debug_print_specific <- function(env = parent.frame(), var_names = c("gt_matrix","gt_std","covariates_matrix","hat")) {
  for (nm in var_names) {
    if (exists(nm, envir = env)) {
      obj <- get(nm, envir = env)
      message("---- ", nm, " ----")
      message(" class: ", paste(class(obj)[1], collapse = ", "))
      if (!is.null(dim(obj))) message(" dims: ", paste(dim(obj), collapse = " x ")) else message(" length: ", length(obj))
      # pretty size
      if (requireNamespace("lobstr", quietly = TRUE)) {
        message(" lobstr::obj_size: ", format(as.numeric(lobstr::obj_size(obj)), big.mark = ","))
      } else {
        message(" object.size: ", format(object.size(obj), units = "auto"))
      }
    } else {
      message("---- ", nm, " not found in env")
    }
  }
  invisible(TRUE)
}
