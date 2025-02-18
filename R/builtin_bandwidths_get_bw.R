#' Get bandwidth functions from string.
#'
#' @keywords internal
#' @param bw_str a string specifying the density of interest.
#' @return a bandwidth function.
get_bw <- function(bw_str) {
  assertthat::assert_that(is.character(bw_str))

  bw <- bw_environment[[bw_str]]

  msg <- paste0("The supplied bandwidth function ('", bw_str, "') is not implemented.")
  assertthat::assert_that(!is.null(bw), msg = msg)

  bw
}


#' Add a new bw to `bw_environment`.
#'
#' @keywords internal
#' @param bw_str A string giving the name of the density.
#' @param bw The bw function.
#' @return None.

add_bw <- function(bw_str, bw) {
  assertthat::assert_that(is.character(bw_str))
  assertthat::assert_that(all(bw_str == make.names(bw_str)),
    msg = "The name of the  bw is not valid. Use a short, valid name. (E.g. kdensity(x, bw = nrd0), where 'nrd0' is a predefined bw function.)"
  )

  func_msg <- paste0("The bw ('", bw_str, "') must be a function.")
  form_msg <- paste0("The bw ('", bw_str, "') must take the arguments 'x', 'kernel', 'start', 'support'.")
  assertthat::assert_that(is.function(bw), msg = func_msg)
  assertthat::assert_that(all(names(formals(bw)) == c("x", "kernel", "start", "support")),
    msg = form_msg
  )

  bw_environment[[bw_str]] <- bw
}



#' Get a bandwidth string when 'bw' is unspecified.
#'
#' @keywords internal
#' @param kernel_str a kernel string
#' @param start_str a parametric start string.
#' @param support the support.
#' @return a bandwidth string.

get_standard_bw <- function(kernel_str, start_str, support) {
  if (kernel_str == "gcopula" & (start_str == "constant" |
    start_str == "uniform")) {
    bw <- "JH"
  } else if (start_str != "constant" & start_str != "uniform") {
    if (!is.null(get_kernel(kernel_str)$sd)) {
      bw <- "RHE"
    } else {
      bw <- "ucv"
    }
  } else if (!is.null(get_kernel(kernel_str)$sd)) {
    bw <- "nrd0"
  } else {
    bw <- "ucv"
  }

  bw
}
