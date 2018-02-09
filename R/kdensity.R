#' Kernel density estimator with parametric start
#'
#' @export
#' @param x the supplied data.
#' @param bw a bandwidth function (NOT IMPLEMENTED) or a double
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
#' @details fill in

kdensity = function(x, adjust = 1, support = NULL, na.rm = TRUE, normalized = TRUE,
                    bw     = c("nrd0",
                               "nrd",
                               "bcv",
                               "ucv",
                               "SJ"),
                    kernel = c("gaussian",
                               "epanechnikov",
                               "rectangular",
                               "triangular",
                               "biweight",
                               "cosine",
                               "optcosine",
                               "tricube",
                               "triweight",
                               "laplace"),
                    start =  c("uniform",
                               "normal",
                               "gamma",
                               "exponential",
                               "inverse_gaussian",
                               "lognormal",
                               "beta",
                               "laplace")) {

  ## We handle the kernel functions. This works by matching the strings and
  ## obtaining the kernel function from outside the 'kdensity' function itself.
  if(!is.list(kernel)) {
    kernel = match.arg(kernel)
    kernel_str = kernel
    kernel = get_kernel(kernel, support)
  } else {
    kernel_str = deparse(substitute(kernel))
  }


  ## Now we handle the parametric start itself. If the parametric start is a
  ## string, we will match the string with pre-supplied functions.

  if(!is.list(start)) {
    start = match.arg(start)
    start_str = start
    start = get_start(start, support)
  } else {
    start_str = deparse(substitute(start))
  }

  ## The support is automatically handled if unspecified.

  if(is.null(support)) {
    support = get_support(start_str, kernel_str)
  }

  ## We continue by checking if the supplied values make sense.
  if(!all(x <= support[2] & x >= support[1])) {
    stop("The supplied data x is not contained in the support: (",
         support[1], ", ", support[2], ").")
  }



  if(!is.element("estimator", names(start)) | !is.element("density", names(start))) {
    stop("The argument 'start' should be either 1.) a string specifying an implemented
          start density or, 2.) a list containing two functions 'density' and 'estimator")
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

  if(!is.numeric(bw)) {
    if(is.character("bw")) {
      bw_str = bw
      bw     = get_bw(bw)(data, kernel, start, support)
    } else {
      bw_str = deparse(substitute(bw))
      bw     = bw(data, kernel, start, support)
    }
  } else {
    bw_str = bw
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
        mean(1/h*kernel$kernel((y-x)/h)*parametric_start_vector(y)/parametric_start_data)
      } else {
        sapply(y, function(y) mean(1/h*kernel$kernel((y-x)/h)*parametric_start_vector(y)/parametric_start_data))
      }
    }

    normalization = tryCatch(integrate(pre_function, lower = support[1], upper = support[2])$value,
                          error = function(e) {
                            stop("Normalization error! The function won't integrate. Try different
                               support?")
                          })

  }


  ## This is the main part of the returned object.

  return_function = function(y) {
    if(length(y) == 1) {
      mean(1/h*kernel$kernel((y-x)/h)*parametric_start_vector(y)/parametric_start_data)/normalization
    } else {
      sapply(y, function(y) {
        parametric_start_vector_y = parametric_start_vector(y)
        1/h*mean(kernel$kernel((y-x)/h)*parametric_start_vector_y/parametric_start_data)/normalization
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

  return_function
}

