#' Get densities and estimators from strings.
#'
#' @keywords internal
#' @param start_str A string specifying the density of interest.
#' @return A list of two functions.

get_start <- function(start_str) {
  assert_(is.character(start_str))

  parametric_start <- starts_environment[[start_str]]

  msg <- paste0("The supplied parametric start ('", start_str, "') is not implemented.")
  assert_(!is.null(parametric_start), msg = msg)

  parametric_start
}

#' Add a new parametric start to `starts_environment`.
#'
#' @keywords internal
#' @param start_str A string giving the name of the density.
#' @param start The parametric start function.
#' @return None.

add_start <- function(start_str, start) {
  assert_(is.character(start_str))
  assert_(all(start_str == make.names(start_str)),
    msg = "The name of the parametric start is not valid. Use a short, valid name. (E.g. kdensity(x, start = gaussian), where gaussian is a predefined start function.)"
  )

  list_msg <- paste0("The parametric start ('", start_str, "') must be a list.")
  assert_(is.list(start), msg = list_msg)

  ## Checks for the right elements in start.
  density_msg <- paste0("The parametric start ('", start_str, "') must contain a function named 'density'.")
  estimator_msg <- paste0("The parametric start ('", start_str, "') must contain a function named 'estimator'.")
  support_msg <- paste0("The parametric start ('", start_str, "') must contain a vector named 'support'.")

  assert_(!is.null(start$density), msg = density_msg)
  assert_(!is.null(start$estimator), msg = estimator_msg)
  assert_(!is.null(start$support), msg = support_msg)

  assign(start_str, start, envir = starts_environment)
}
