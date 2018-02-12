#' Kernel density estimator with parametric start
#'
#' \code{kdensity} computes a kernel density for univariate data. It supports
#' asymmetric kernels and parametric starts through the \code{kernel} and
#' \code{start} arguments.
#'
#' @export
#'
#' @param x Numeric vector; the data.
#'
#' @param bw A bandwidth function. Can be either a string, a custom made
#' function, or a double.
#'
#' @param adjust An adjustment constant, so that \code{h = adjust*bw*sd}, where \code{sd}
#' varies acccording to the chosen kernel.
#'
#' @param kernel The kernel function. Can be a string or a custom made list
#' containg a function \code{kernel}, the standard deviation of the kernel,
#' \code{sd}, and domain of definition for the kernel, \code{support}.
#'
#' @param start Choice of parametric start. Can be a string or be supplied via a l
#' ist containing a \code{density} function, an \code{estimator} function,
#' and a \code{support} tuple.
#'
#' @param support The support of the data. Must be compatible with the supplied
#' \code{x} and the supplied \code{start} and \code{kernel}. Is used to find the
#' normalization constant, see \code{normalized}.
#'
#' @param normalized Logical; if \code{TRUE}, the density is normalized.
#'
#' @param na.rm Logical; if \code{TRUE}, \code{NA}s will be removed from \code{x}.
#'
#' @return \code{kdensity} Returns a function object of \code{\link[base]{class}} "kdensity".
#'
#' @details If \code{normalized} is \code{FALSE} and \code{start != "uniform"}, the resulting
#' density will not integrate to 1 in general.
#'
#'   \strong{Bandwidth functions}: Bandwidth functions can either be specified by a string,
#'   be user made, or be a fixed double. If the argument is a string, it must fully match
#'   one of the implemented bandwidt functions. From the package \code{stats}, the bandwidth
#'   functions are \code{nrd0}, \code{nrd}, \code{bcv}, \code{ucv}, and \code{SJ} are avaiable,
#'   see \code{\link[stats]{bandwidth}}. They intended for use with a 'uniform' start and the 'gaussian' kernel, but
#'   work well for the other symmetric kernels as well. Implemented bandwidth functions for asymmetric
#'   kernels are: \code{JH} for the Gaussian copula estimator. For parametric starts, \code{RHE} is an
#'   Hermite expansion reference rule. Bandwidth functions for asymmetric kernels and parametric
#'   starts are documented in \code{\link{bandwidth_functions}}.
#'
#'   \strong{Kernel functions}: Kernel functions can either be specified by a string or
#'   be user made. If the argument is a string, it must fully match one of the implemented
#'   kernels. Available symmetric kernels are \code{gaussian} (or \code{normal}), \code{epanechnikov},
#'   \code{rectangular} (or \code{uniform}), \code{triangular}, \code{biweight}, \code{triweight},
#'   \code{tricube}, \code{cosine}, \code{optcosine}, and \code{laplace}.
#'   See \code{\link[stats]{density}} for more details.
#'
#'  The implemented asymmetric kernels are:
#'   \itemize{
#'     \item \code{gcopula}. The Gaussian copula KDE, used for data on the unit
#'      interval. Described in Jones & Henderswon.
#'     \item \code{gamma} and \code{gamma_biased}. Gamma kernels for data on the
#'     positive half-line, \code{c(0, Inf)}. They are described in Chen.
#'   }
#'
#'   \strong{Parametric starts}: The following parametric starts are supported:
#'   \code{uniform} (or \code{constant}), \code{normal}, \code{gamma},
#'   \code{exponential}, \code{inverse_gaussian}, \code{lognormal}, \code{beta},
#'   and \code{laplace}. Their parameters are estimated by maximum likelihood.
#'   The default value is \code{uniform}, which corresponds to ordinary kernel
#'   density estimation.
#' @seealso The \code{stats} package function \code{\link[stats]{density}}. For
#'   bandwidth selection documentation, see \code{\link{bandwidth_seletor}}.
#' @references Hjort, Nils Lid, and Ingrid K. Glad. "Nonparametric density estimation with a parametric start." The Annals of Statistics (1995): 882-904.
#'
#'   Jones, M. C., and D. A. Henderson. "Miscellanea kernel-type density estimation on the unit interval." Biometrika 94.4 (2007): 977-984.
#'
#'   Chen, Song Xi. "Probability density function estimation using gamma kernels." Annals of the Institute of Statistical Mathematics 52.3 (2000): 471-480.
#'
#'   Silverman, Bernard W. Density estimation for statistics and data analysis. Vol. 26. CRC press, 1986.
#'
#'  @examples
#' ## Use gamma kernels to model positive data, the concentration of
#' ## theophylline
#'
#' concentration = Theoph$conc + 0.001
#' plot(kdensity(concentration, start = "gamma", kernel = "gamma", adjust = 1/3),
#'      ylim = c(0, 0.15), lwd = 2, main = "Concentration of theophylline")
#' lines(kdensity(concentration, start = "gamma", kernel = "gaussian"),
#'       lty = 2, col = "grey", lwd = 2)
#' lines(kdensity(concentration, start = "gaussian", kernel = "gaussian"),
#'       lty = 3, col = "blue", lwd = 2)
#' lines(kdensity(concentration, start = "gaussian", kernel = "gamma", adjust = 1/3),
#'       lty = 4, col = "red", lwd = 2)
#' rug(concentration)
#'
#' ## Using a density and and estimator from another package.
#'
#' skew_hyperbolic = list(
#'   density   = SkewHyperbolic::dskewhyp,
#'   estimator = function(x) SkewHyperbolic::skewhypFit(dat,printOut = FALSE)$param,
#'   support   = c(-Inf, Inf)
#' )
#'
#' kde = kdensity(diff(LakeHuron), start = skew_hyperbolic)
#' plot(kde, lwd = 2, col = "blue",
#'      main = "Annual differences in water level (ft) of Lake Huron, 1875 - 1972")
#' lines(kde, plot_start = TRUE, lty = 2, lwd = 2) # Plots the skew hyperbolic density.
#' rug(diff(LakeHuron))

kdensity = function(x, bw = NULL, adjust = 1, kernel = NULL, start = NULL,
                    support = NULL, na.rm = FALSE, normalized = TRUE)
 {

  data.name = deparse(substitute(x))
  has.na = any(is.na(x))

  if(has.na) {
    if(!na.rm) stop("x contains NAs and na.rm = FALSE.")
    x = x[!is.na(x)]
  }

  ## The case of bw == Inf is special! In this case, we return the parametric
  ## start itself.

  if(!is.null(bw)) {
    if(is.numeric(bw)) {
      if(bw == Inf) {

        if(!is.list(start)) {
          msg = "bw = Inf does not work with a uniform start."
          assertthat::assert_that(!is.null(start), start != "uniform", msg = msg)
        }

        kss_list = get_kernel_start_support(NULL, start, NULL)
        start_str = ifelse(!is.list(start), kss_list$start_str, deparse(substitute(start)))
        start = kss_list$start

        parameters = start$estimator(x)
        parametric_start = start$density

        return_function = function(y) {
          sapply(y, function(y) {
            do.call(parametric_start, as.list(c("x" = y, parameters)))
          })
        }

        class(return_function) = "kdensity"
        attr(return_function, "bw_str")    = Inf
        attr(return_function, "bw")        = Inf
        attr(return_function, "kernel")    = "none"
        attr(return_function, "start")     = start_str
        attr(return_function, "support")   = start$support
        attr(return_function, "adjust")    = 1
        attr(return_function, "n")         = length(x)
        attr(return_function, "h")         = Inf
        attr(return_function, "data.name") = deparse(substitute(x))
        attr(return_function, "has.na")    = any(is.na(x))
        attr(return_function, "call")      = match.call()
        attr(return_function, "range")     = c(min(x), max(x))
        attr(return_function, "estimates") = parameters
        return(return_function)
      }
    }
  }


  ## Now we massage and handle the combinations of kernel, start and support.
  kss_list = get_kernel_start_support(kernel, start, support)

  start_str = ifelse(!is.list(start), kss_list$start_str, deparse(substitute(start)))
  kernel_str = ifelse(!is.list(kernel), kss_list$kernel_str, deparse(substitute(kernel)))

  kernel = kss_list$kernel
  start = kss_list$start
  support = kss_list$support

  ## Tests for incompabibilities in the supplied values.
  support_compatible(kernel, start, support)

  ## We continue by checking if the supplied values make sense.
  if(!all(x <= support[2] & x >= support[1])) {
    stop("The supplied data x is not contained in the support: (",
         support[1], ", ", support[2], ").")
  }

  parameters = start$estimator(x)
  parametric_start = start$density
  # Name of the variable where the density is evaluated. Typically x.
  x_name = names(formals(start$density))[1]

  parametric_start_vector = function(data) {
    sapply(data, function(datum) {
      arguments = as.list(c("x" = datum, parameters))
      names(arguments)[1] = x_name
      do.call(parametric_start, arguments)
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
  sd = ifelse(!is.null(kss_list$kernel$sd), kss_list$kernel$sd, 1)
  h = bw*adjust*sd

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
        sapply(y, function(y) mean(1/h*kernel$kernel(y, x, h)/parametric_start_data)*parametric_start_vector(y))
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
  attr(return_function, "data.name") = data.name
  attr(return_function, "has.na")    = has.na
  attr(return_function, "call")      = match.call()
  attr(return_function, "range")     = c(min(x), max(x))
  attr(return_function, "estimates") = parameters

  return_function
}

