checkpointed_window_source_file <- tryCatch(
  normalizePath(sys.frame(1L)$ofile, mustWork = TRUE),
  error = function(condition) {
    file_argument <- grep(
      "^--file=",
      commandArgs(trailingOnly = FALSE),
      value = TRUE
    )
    if (length(file_argument) != 1L) {
      return(normalizePath(
        "R/scripts/run_checkpointed_window_susie.R",
        mustWork = TRUE
      ))
    }
    script_path <- sub("^--file=", "", file_argument)
    script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
    normalizePath(script_path, mustWork = TRUE)
  }
)

FunctionPath <- Sys.getenv("SUSIER_FUNCTIONS_PATH", unset = "/opt/r/lib")
if (!file.exists(file.path(FunctionPath, "CheckpointedWindowSusieFunctions.R"))) {
  FunctionPath <- file.path(
    dirname(dirname(checkpointed_window_source_file)),
    "utils"
  )
}

source(file.path(FunctionPath, "InitFunctions.R"))
source(file.path(FunctionPath, "SusieFunctions.R"))
source(file.path(FunctionPath, "CheckpointedWindowSusieFunctions.R"))
source(file.path(FunctionPath, "CheckpointStore.R"))

checkpointed_window_option_list <- function() {
  list(
    optparse::make_option(
      "--window-id",
      dest = "window_id",
      type = "character",
      help = "Prepared window identifier."
    ),
    optparse::make_option(
      "--window-dosage",
      dest = "window_dosage",
      type = "character",
      help = "Prepared window dosage TSV."
    ),
    optparse::make_option(
      "--window-phenotypes",
      dest = "window_phenotypes",
      type = "character",
      help = "Ordered-window phenotype manifest TSV."
    ),
    optparse::make_option(
      "--phenotype-data",
      dest = "phenotype_data",
      type = "character",
      help = "Prepared phenotype data file."
    ),
    optparse::make_option(
      "--covariate-files",
      dest = "covariate_files",
      type = "character",
      help = "Comma-separated covariate files."
    ),
    optparse::make_option(
      "--covariate-modalities",
      dest = "covariate_modalities",
      type = "character",
      help = "Comma-separated covariate modality labels."
    ),
    optparse::make_option(
      "--keep-samples",
      dest = "keep_samples",
      type = "character",
      default = NULL,
      help = "Optional sample allowlist."
    ),
    optparse::make_option(
      "--checkpoint-root",
      dest = "checkpoint_root",
      type = "character",
      help = "Filesystem path or GCS URI for checkpoints."
    ),
    optparse::make_option(
      "--output-dir",
      dest = "output_dir",
      type = "character",
      help = "Local directory for completed window outputs."
    )
  )
}

checkpointed_window_option_parser <- function() {
  optparse::OptionParser(
    usage = "%prog [options]",
    description = "Run ordered checkpointed univariate SuSiE for one window.",
    option_list = checkpointed_window_option_list()
  )
}

checkpointed_required_option <- function(value, flag) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(trimws(value))) {
    stop("Required option is absent: ", flag, call. = FALSE)
  }
  value
}

parse_checkpointed_comma_array <- function(value, flag) {
  value <- checkpointed_required_option(value, flag)
  if (grepl("(^,|,$|,,)", value)) {
    stop(flag, " cannot contain empty members.", call. = FALSE)
  }
  members <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  if (length(members) == 0L || any(!nzchar(members))) {
    stop(flag, " cannot contain empty members.", call. = FALSE)
  }
  members
}

parse_checkpointed_window_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  options <- optparse::parse_args(
    checkpointed_window_option_parser(),
    args = args,
    positional_arguments = FALSE
  )
  covariate_files <- parse_checkpointed_comma_array(
    options$covariate_files,
    "--covariate-files"
  )
  covariate_modalities <- parse_checkpointed_comma_array(
    options$covariate_modalities,
    "--covariate-modalities"
  )
  if (length(covariate_files) == 0L ||
      length(covariate_files) != length(covariate_modalities)) {
    stop(
      "--covariate-files and --covariate-modalities must have equal nonzero lengths.",
      call. = FALSE
    )
  }

  list(
    window_id = checkpointed_required_option(options$window_id, "--window-id"),
    window_dosage = checkpointed_required_option(
      options$window_dosage,
      "--window-dosage"
    ),
    window_phenotypes = checkpointed_required_option(
      options$window_phenotypes,
      "--window-phenotypes"
    ),
    phenotype_data = checkpointed_required_option(
      options$phenotype_data,
      "--phenotype-data"
    ),
    covariate_files = covariate_files,
    covariate_modalities = covariate_modalities,
    keep_samples = if (is.null(options$keep_samples)) {
      NULL
    } else {
      checkpointed_required_option(options$keep_samples, "--keep-samples")
    },
    checkpoint_root = checkpointed_required_option(
      options$checkpoint_root,
      "--checkpoint-root"
    ),
    output_dir = checkpointed_required_option(options$output_dir, "--output-dir"),
    wrapper_path = checkpointed_window_source_file
  )
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  config <- parse_checkpointed_window_cli(args)
  run_checkpointed_window(config)
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  main()
}
