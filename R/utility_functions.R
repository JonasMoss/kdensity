assert_ <- function(..., msg = "Assertion failed.") {
  conditions <- list(...)
  ok <- all(vapply(conditions, function(condition) isTRUE(all(condition)), logical(1)))

  if (!ok) {
    stop(msg, call. = FALSE)
  }

  invisible(TRUE)
}
