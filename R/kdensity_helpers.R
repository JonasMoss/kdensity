#' Fill in missing kernel, start or support given the supplied values.
#'
#' This function takes the supplied values of kernel, start, and support
#' and fills in the non-supplied ones. It also handles inconsistencies,
#' such as providing a support on (-Inf, Inf) but a kernel on (0, Inf).
#'
#' The `kernel` and `start` parameters are either strings or
#' adhering to the kernel/start list structure. `support` is a
#' numeric vector of length two.
#'
#' @keywords internal
#' @param kernel Supplied kernel; string or list.
#' @param start Supplied parametric start; string or list.
#' @param support Binary vector.
#' @return a list with members kernel, kernel_str, start, start_str,
#' and support.

get_kernel_start_support <- function(kernel, start, support) {
  ## First we handle the parametric start. This is easy, since its default
  ## value equals uniform.

  if (!is.null(start)) {
    start_str <- start
    start <- get_start(start)
  } else {
    start_str <- "uniform"
    start <- get_start("uniform")
  }

  if (is.null(kernel) & is.null(support)) {
    ## No arguments will give you the stats::density behaviour by default.
    return(list(
      kernel = get_kernel("gaussian"),
      kernel_str = "gaussian",
      start = start,
      start_str = start_str,
      support = start$support
    ))
  }


  ## What happens when kernel is also filled in?

  if (is.null(kernel) & !is.null(support)) {
    ## I match a kernel with a support.
    if (all(support == c(-Inf, Inf))) {
      kernel <- "gaussian"
    } else if (all(support == c(0, 1))) {
      kernel <- "gcopula"
    } else if (all(support == c(0, Inf))) {
      kernel <- "gamma"
    } else {
      message("The supplied support does not exactly match any standard kernels.")
      l <- support[1]
      u <- support[2]
      if (l < 0) {
        kernel <- "gaussian"
      } else {
        if (u < 1) {
          kernel <- "gcopula"
        } else {
          kernel <- "gamma"
        }
      }
    }
  }

  ## The next step is to fill in the kernel if it is non-NULL:
  if (!is.null(kernel)) {
    kernel_str <- kernel
    kernel <- get_kernel(kernel)
  }

  if (!is.null(kernel) & is.null(support)) {
    ## The resulting support must be compatible with both the kernel
    ## and the start.
    l <- max(start$support[1], kernel$support[1])
    u <- min(start$support[2], kernel$support[2])
    return(list(
      kernel = kernel,
      kernel_str = kernel_str,
      start = start,
      start_str = start_str,
      support = c(l, u)
    ))
  }

  return(list(
    kernel = kernel,
    kernel_str = kernel_str,
    start = start,
    start_str = start_str,
    support = support
  ))
}


#' Checks compatibility between supports.
#'
#' The supplied support must never be larger than the support of
#' the parametric start / kernel.
#' @keywords internal
#' @param kernel,start,support The kernel, start and support to check.
#' @return None.
support_compatible <- function(kernel, start, support) {
  assertthat::assert_that(kernel$support[1] <= support[1],
    msg =
      "The lower end point of the support is smaller than the lower end point of the 'kernel support'."
  )

  assertthat::assert_that(kernel$support[2] >= support[2],
    msg =
      "The upper end point of the support is larger than the upper end point of the 'kernel support'."
  )

  assertthat::assert_that(start$support[1] <= support[1],
    msg =
      "The lower end point of the support is smaller than the lower end point of the 'start support'."
  )

  assertthat::assert_that(start$support[2] >= support[2],
    msg =
      "The upper end point of the support is larger than the upper end point of the 'start support'."
  )
}
