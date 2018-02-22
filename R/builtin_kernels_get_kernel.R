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
