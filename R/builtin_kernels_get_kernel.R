#' Helper function that gets a kernel function for kdensity.
#'
#' @keywords internal
#' @param kernel_str a string specifying which kernel to use.
#' @return a kernel function of the format k(u) with integral normalized
#' to 1.

get_kernel <- function(kernel_str) {
  assertthat::assert_that(is.character(kernel_str))

  kernel <- kernel_environment[[kernel_str]]

  msg <- paste0("The supplied kernel ('", kernel_str, "') is not implemented.")
  assertthat::assert_that(!is.null(kernel), msg = msg)

  kernel
}

#' Add a new kernel to `kernels_environment`.
#'
#' @keywords internal
#' @param kernel_str A string giving the name of the density.
#' @param kernel The kernel function.
#' @return None.

add_kernel <- function(kernel_str, kernel) {
  assertthat::assert_that(is.character(kernel_str))
  assertthat::assert_that(all(kernel_str == make.names(kernel_str)),
    msg = "The name of the kernel is not valid. Use a short, valid name. (E.g. kdensity(x, kernel = gaussian), where gaussian is a predefined kernel function.)"
  )

  list_msg <- paste0("The kernel ('", kernel_str, "') must be a list.")
  assertthat::assert_that(is.list(kernel), msg = list_msg)

  ## Checks for the right elements in kernel.
  density_msg <- paste0("The kernel ('", kernel_str, "') must contain a function named 'kernel'.")
  support_msg <- paste0("The kernel ('", kernel_str, "') must contain a function named 'support'.")
  assertthat::assert_that(!is.null(kernel$kernel), msg = density_msg)
  assertthat::assert_that(!is.null(kernel$support), msg = support_msg)

  kernel_environment[[kernel_str]] <- kernel
}
