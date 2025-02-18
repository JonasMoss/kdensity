#' Get densities and estimators from strings.
#'
#' @keywords internal
#' @param start_str A string specifying the density of interest.
#' @return A list of two functions.

get_start <- function(start_str) {
  assertthat::assert_that(is.character(start_str))

  parametric_start <- starts_environment[[start_str]]

  msg <- paste0("The supplied parametric start ('", start_str, "') is not implemented.")
  assertthat::assert_that(!is.null(parametric_start), msg = msg)

  parametric_start
}

#' Add a new parametric start to `starts_environment`.
#'
#' @keywords internal
#' @param start_str A string giving the name of the density.
#' @param start The parametric start function.
#' @return None.

add_start <- function(start_str, start) {
  assertthat::assert_that(is.character(start_str))
  assertthat::assert_that(all(start_str == make.names(start_str)),
    msg = "The name of the parametric start is not valid. Use a short, valid name. (E.g. kdensity(x, start = gaussian), where gaussian is a predefined start function.)"
  )

  list_msg <- paste0("The parametric start ('", start_str, "') must be a list.")
  assertthat::assert_that(is.list(start), msg = list_msg)

  ## Checks for the right elements in start.
  density_msg <- paste0("The parametric start ('", start_str, "') must contain a function named 'density'.")
  estimator_msg <- paste0("The parametric start ('", start_str, "') must contain a function named 'estimator'.")
  support_msg <- paste0("The parametric start ('", start_str, "') must contain a vector named 'support'.")

  assertthat::assert_that(!is.null(start$density), msg = density_msg)
  assertthat::assert_that(!is.null(start$estimator), msg = estimator_msg)
  assertthat::assert_that(!is.null(start$support), msg = support_msg)

  assign(start_str, start, envir = starts_environment)
}
