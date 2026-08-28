expect_error_message <- function(code, expected_message) {
  error_message <- tryCatch(
    {
      force(code)
      NULL
    },
    error = function(error) conditionMessage(error)
  )
  if (is.null(error_message)) {
    stop("Expected an error containing: ", expected_message, call. = FALSE)
  }
  if (!grepl(expected_message, error_message, fixed = TRUE)) {
    stop(
      "Expected error containing '", expected_message,
      "' but received: ", error_message,
      call. = FALSE
    )
  }
  invisible(TRUE)
}

expect_identical_value <- function(actual, expected, label = "value") {
  if (!identical(actual, expected)) {
    stop(
      "Unexpected ", label, ".\nExpected: ", paste(capture.output(dput(expected)), collapse = ""),
      "\nActual: ", paste(capture.output(dput(actual)), collapse = ""),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

run_named_test <- function(name, code) {
  force(code)
  cat("PASS ", name, "\n", sep = "")
  invisible(TRUE)
}
