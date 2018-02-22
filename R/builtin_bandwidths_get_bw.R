#' Get bandwidth functions from string.
#'
#' @param bw_str a string specifying the density of interest.
#' @return a bandwidth function.
get_bw = function(bw_str) {

  assertthat::assert_that(is.character(bw_str))

  bw = bw_environment[[bw_str]]

  if(is.null(bw)) {
    if(exists(bw_str)) {
      parametric_start = get(bw_str)
    } else {
      stop(paste0("The supplied kernel ('", bw_str,"') is not implemented."))
    }
  }

  bw

}


#' Add a new bw to \code{bw_environment}.
#'
#' @param bw_str A string giving the name of the density.
#' @param bw The bw function.
#' @return None.

add_bw = function(bw_str, bw) {
  assertthat::assert_that(is.character(bw_str))
  assertthat::assert_that(all(bw_str == make.names(bw_str)),
                          msg = "The name of the  bw is not valid. Use a short, valid name. (E.g. kdensity(x, bw = nrd0), where 'nrd0' is a predefined bw function.)")

  func_msg = paste0("The bw ('", bw_str, "') must be a function.")
  form_msg = paste0("The bw ('", bw_str, "') must take the arguments 'x', 'kernel', 'start', 'support'.")
  assertthat::assert_that(is.function(bw), msg = func_msg)
  assertthat::assert_that(all(formals(bw) == c("x", "kernel", "start", "support")),
                          msg = form_msg)

  bw_environment[[bw_str]] = bw
}



#' Get a bandwidth string when 'bw' is unspecified.
#'
#' @param kernel_str a kernel string
#' @param start_str a parametric start string.
#' @param support the support.
#' @return a bandwidth string.

get_standard_bw = function (kernel_str, start_str, support) {
  if(start_str != "constant" & start_str != "uniform") {
    if(kernel_str != "gaussian" & kernel_str != "normal") {
      if(kernel_str == "gcopula") {
        bw = "JH"
      } else {
        bw = "ucv"
      }
    } else {
      bw = "RHE"
    }
  } else {
    bw ="nrd0"
  }
  bw
}
