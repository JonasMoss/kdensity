assert_ <- function(..., msg = "Assertion failed.") {
  conditions <- list(...)
  ok <- all(vapply(conditions, function(condition) isTRUE(all(condition)), logical(1)))

  if (!ok) {
    stop(msg, call. = FALSE)
  }

  invisible(TRUE)
}

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
  assert_(kernel$support[1] <= support[1],
    msg =
      "The lower end point of the support is smaller than the lower end point of the 'kernel support'."
  )

  assert_(kernel$support[2] >= support[2],
    msg =
      "The upper end point of the support is larger than the upper end point of the 'kernel support'."
  )

  assert_(start$support[1] <= support[1],
    msg =
      "The lower end point of the support is smaller than the lower end point of the 'start support'."
  )

  assert_(start$support[2] >= support[2],
    msg =
      "The upper end point of the support is larger than the upper end point of the 'start support'."
  )
}

## Resolve the `bw` argument of kdensity() into a numeric bandwidth and a
## reporting string. `bw` may be NULL (use the standard selector), a string
## naming a built-in selector, a custom function, or a numeric value.
## `bw_label` is `deparse(substitute(bw))` captured by the caller, used only to
## label a custom bandwidth function.
resolve_bandwidth <- function(bw, bw_label, x, kernel_str, start_str, support) {
  if (is.null(bw)) {
    bw <- get_standard_bw(kernel_str, start_str, support)
  }

  if (!is.numeric(bw)) {
    if (is.character(bw)) {
      bw_str <- bw
      bw <- get_bw(bw)(x, kernel_str, start_str, support)
    } else {
      bw_str <- bw_label
      bw <- bw(x, kernel_str, start_str, support)
    }
  } else {
    bw_str <- "user supplied"
    if (bw == Inf) {
      msg <- "bw = Inf does not work with a uniform start."
      assert_(start_str != "uniform", start_str != "constant", msg = msg)
    }
  }

  list(bw = bw, bw_str = bw_str)
}

## Compute the normalization constant for the unnormalized density. Uses a
## known closed-form normalizer when available, otherwise numerical integration
## with a relative-error tolerance check.
resolve_normalization <- function(unnormalized_density, start, kernel, support,
                                  h, x, estimates, tolerance) {
  normalization <- known_normalization(
    start = start,
    kernel = kernel,
    support = support,
    h = h,
    x = x,
    estimates = estimates
  )

  if (!is.null(normalization)) {
    return(normalization)
  }

  integral <- tryCatch(
    stats::integrate(unnormalized_density,
      lower = support[1],
      upper = support[2]
    ),
    error = function(e) {
      stop(paste0(
        "Normalization error: The function will not integrate.",
        "Two common causes are: 1.) The kernel is non-smooth, ",
        "try a smooth kernel if possible. 2.) The supplied ",
        "support is incorrect."
      ))
    }
  )

  msg <- "The normalization constant has too large relative error.
   1.) try to normalize the data, or
   2.) change kernel, or
   3.) increase the 'tolerance' parameter for kdensity"

  assert_(!is.nan(integral$abs.error / integral$value),
    msg = msg
  )
  assert_(integral$abs.error / integral$value < tolerance,
    msg = msg
  )

  integral$value
}
