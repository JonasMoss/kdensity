#' Helper function that gets a kernel function for kdensity.
#'
#' @param kernel_str a string specifying which kernel to use.
#' @return a kernel function of the format k(u) with integral normalized
#' to 1.

get_kernel = function(kernel_str) {

  assertthat::assert_that(is.character(kernel_str))

  kernel = kernel_environment[[kernel_str]]

  if(is.null(kernel)) {
    if(exists(kernel_str)) {
      parametric_start = get(kernel_str)
    } else {
      stop(paste0("The supplied kernel ('", kernel_str,"') is not implemented."))
    }
  }

  kernel
}


#' Add a new parametric start to \code{starts_environment}.
#'
#' @param start_str A string giving the name of the density.
#' @param start The parametric start function.
#' @return None.

add_start = function(start_str, start) {
  assertthat::assert_that(is.character(start_str))
  assertthat::assert_that(all(start_str == make.names(start_str)),
                          msg = "The name of the parametric start is not valid. Use a short, valid name. (E.g. kdensity(x, start = gaussian), where gaussian is a predefined start function.)")

  list_msg = paste0("The parametric start ('", start_str, "') must be a list.")
  assertthat::assert_that(is.list(start), msg = list_msg)

  ## Checks for the right elements in start.
  density_msg = paste0("The parametric start ('", start_str, "') must contain a function named 'density'.")
  estimator_msg = paste0("The parametric start ('", start_str, "') must contain a function named 'estimator.")
  support_msg = paste0("The parametric start ('", start_str, "') must contain a function named 'support'.")

  assertthat::assert_that(!is.null(start$density), msg = density_msg)
  assertthat::assert_that(!is.null(start$estimator), msg = estimator_msg)
  assertthat::assert_that(!is.null(start$support), msg = support_msg)

  starts_environment[[start_str]] = start
}


#' Add a new kernel to \code{kernels_environment}.
#'
#' @param kernel_str A string giving the name of the density.
#' @param kernel The kernel function.
#' @return None.

add_kernel = function(kernel_str, kernel) {
  assertthat::assert_that(is.character(kernel_str))
  assertthat::assert_that(all(kernel_str == make.names(kernel_str)),
                          msg = "The name of the parametric kernel is not valid. Use a short, valid name. (E.g. kdensity(x, kernel = gaussian), where gaussian is a predefined kernel function.)")

  list_msg = paste0("The kernel ('", kernel_str, "') must be a list.")
  assertthat::assert_that(is.list(kernel), msg = list_msg)

  ## Checks for the right elements in kernel.
  density_msg = paste0("The kernel ('", kernel_str, "') must contain a function named 'kernel'.")
  support_msg = paste0("The kernel ('", kernel_str, "') must contain a function named 'support'.")
  assertthat::assert_that(!is.null(kernel$kernel), msg = density_msg)
  assertthat::assert_that(!is.null(kernel$support), msg = support_msg)

  kernel_environment[[kernel_str]] = kernel
}

