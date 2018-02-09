#' Kernel density estimator with parametric start
#'
#' @export
#' @param x the supplied data.
#' @param bw a bandwidth function or a double.
#' @param adjust an adjustment constant, so that h = adjust*bw*m, where m
#' varies acccording to the chosen kernel.
#' @param kernel the kernel function, chosen from a list of alternatives
#' @param start choice of parametric start. Can either be chosen from a vector
#' of strings, or be supplied via a list containing a \code{density} function
#' and an \code{estimator} function.
#' @param support the support of the data. Must be compatible with the supplied
#' \code{x} and the supplied \code{start} and \code{kernel}
#' @param normalized should the density be normalized?
#' @param na.rm should \code{NA} values be removed?
#' @return A function of class kdensity.
#' @details
#'
#' Bandwidth functions: The available bandwidth functions are "nrd0", "nrd",
#' "bcv", "ucv", and "SJ" from the stats package. They intended for use with a
#' 'uniform' start and the 'gaussian' kernel, but work well for the other
#' symmetric kernels as well.
#'
#' The function "JH" is designed for use with the 'gcopula' kernel.
#'
#' The standard bandwidth function is 'nrd0' when the start is uniform and
#' the kernel is symmetric, following stats::density.
#'
#' When the kernel "gcopula" is chosen, the standard bandwidth function is
#' "JH".
#'
#' (NOT IMPLEMENTED) When a symmetric kernel is chosen and start != uniform,
#' the alternatives are:
#'
#'
kdensity = function(x, kernel = NULL, start = NULL, bw = NULL, adjust = 1,
                    support = NULL, na.rm = FALSE, normalized = TRUE)
 {

  ## Now we massage and handle the combinations of kernel, start and support.
  kss_list = get_kernel_start_support(kernel, start, support)

  kernel = kss_list$kernel
  kernel_str = kss_list$kernel_str
  start = kss_list$start
  start_str = kss_list$start_str
  support = kss_list$support

  ## We continue by checking if the supplied values make sense.
  if(!all(x <= support[2] & x >= support[1])) {
    stop("The supplied data x is not contained in the support: (",
         support[1], ", ", support[2], ").")
  }

  parameters = start$estimator(x)
  parametric_start = start$density

  parametric_start_vector = function(data) {
    sapply(data, function(datum) {
      do.call(parametric_start, as.list(c("x" = datum, parameters)))
      })
  }

  ## Takes care of the bandwidth. Can be either a double, a string, or a
  ## function taking the arguments data, kernel, start support.

  if(is.null(bw)) {
    bw = get_standard_bw(kernel_str, start_str, support)
  }

  if(!is.numeric(bw)) {
    if(is.character(bw)) {
      bw_str = bw
      bw     = get_bw(bw)(x, kernel_str, start_str, support)
    } else {
      bw_str = deparse(substitute(bw))
      bw     = bw(data, kernel_str, start_str, support)
    }
  } else {
    bw_str = "user supplied"
  }

  ## The parameter h is computed. The basic bandwidth is h = bw*adjust for the
  ## normal kernel, and is adjusted for all the other kernels so that the sd
  ## of the kernel equals h.
  h = bw*adjust*kernel$sd

  ## The denominator can be computed once and for all.
  parametric_start_data = parametric_start_vector(x)

  ## Now we handle normalization. If the supplied start density is uniform,
  ## there is nothing to do. If it isn't uniform, the R-function integrate is
  ## used to find the normalization constant unless normalized is FALSE.

  normalization = 1

  if(normalized & !(start_str == "uniform" & all(support == c(-Inf, Inf)))) {

    pre_function = function(y) {
      if(length(y) == 1) {
        mean(1/h*kernel$kernel(y, x, h)*parametric_start_vector(y)/parametric_start_data)
      } else {
        sapply(y, function(y) mean(1/h*kernel$kernel(y, x, h)*parametric_start_vector(y)/parametric_start_data))
      }
    }

    normalization = tryCatch(integrate(pre_function, lower = support[1], upper = support[2])$value,
                          error = function(e) {
                            stop("Normalization error: The function will not integrate. Two common causes are: 1.) The kernel is non-smooth, try a smooth kernel if possible. 2.) The supplied support is incorrect.")
                          })

  }


  ## This is the main part of the returned object.

  return_function = function(y) {
    if(length(y) == 1) {
      mean(1/h*kernel$kernel(y, x, h)*parametric_start_vector(y)/parametric_start_data)/normalization
    } else {
      sapply(y, function(y) {
        parametric_start_vector_y = parametric_start_vector(y)
        1/h*mean(kernel$kernel(y, x, h)*parametric_start_vector_y/parametric_start_data)/normalization
      })
    }

  }


  ## Finally, we assign the return_function its class and required attributes.
  class(return_function) = "kdensity"
  attr(return_function, "bw_str")    = bw_str
  attr(return_function, "bw")        = bw
  attr(return_function, "kernel")    = kernel_str
  attr(return_function, "start")     = start_str
  attr(return_function, "support")   = support
  attr(return_function, "adjust")    = adjust
  attr(return_function, "n")         = length(x)
  attr(return_function, "h")         = h
  attr(return_function, "data.name") = deparse(substitute(x))
  attr(return_function, "has.na")    = any(is.na(x))
  attr(return_function, "call")      = match.call()
  attr(return_function, "range")     = c(min(x), max(x))
  attr(return_function, "estimates") = parameters

  return_function
}

