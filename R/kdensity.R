#' Parametrically guided kernel density estimation
#'
#' `kdensity` computes a parametrically guided kernel density estimate
#' for univariate data. It supports asymmetric kernels and parametric starts
#' through the `kernel` and `start` arguments.
#'
#' The default values for `bw`, `kernel`, `start`, and
#' `support` are interdependent, and are chosen to make sense. E.g.,
#' the default value for `support` when `start = beta` is
#' `c(0, 1)`.
#'
#' The `start` argument defaults to `uniform`, which corresponds
#' to ordinary kernel density estimation. The typical default value for
#' `kernel` is `gaussian`.
#'
#' @export
#'
#' @param x Numeric vector containing the data.
#'
#' @param bw A bandwidth function. Can be either a string, a custom-made
#' function, or a double. The supported bandwidth functions are documented
#' in [bandwidths()].
#'
#' @param adjust An adjustment constant, so that `h = adjust*bw*sd`, where `sd`
#' varies with the chosen kernel.
#'
#' @param kernel The kernel function. Can be chosen from the list of built-in
#' kernels or be custom-made. See [kernels()] for details.
#'
#' @param start Parametric start. Can be chosen from the list of built-in
#' parametric starts or be custom-made. See [parametric_starts()] for
#' details.
#'
#' @param support The support of the data. Must be compatible with the supplied
#' `x` and the supplied `start` and `kernel`. Is used to find the
#' normalization constant, see `normalized`.
#'
#' @param na.rm Logical; if `TRUE`, `NA`s will be removed from `x`.
#'
#' @param normalized Logical; if `TRUE`, the density is normalized.
#'
#' @param tolerance Numeric; the relative error to tolerate in normalization.
#'
#' @return `kdensity` returns an S3 function object of
#' [base::class()] "kdensity". This is a callable function with the
#' following elements, accessible by '$':
#' \describe{
#'   \item{`x`}{The data supplied in `x`.}
#'   \item{`bw_str, bw, adjust, h`}{The bandwidth function, the resulting
#'               bandwidth, the `adjust` argument, and the adjusted
#'               bandwidth.}
#'   \item{`kernel_str, kernel, start, start_str, support`}{Name of the kernel,
#'               the kernel object, name of the parametric start, the start object,
#'               and the support of the density.}
#'   \item{`data.name, n, range, has.na, na.rm, normalized`}{Name of the data, number of
#'               observations, the range of the data, whether the data
#'               `x` contained `NA` values, whether na.rm is `TRUE`
#'               or not, and whether the density is normalized.}
#'   \item{`call`}{The `call` to `kdensity`.}
#'   \item{`estimates`}{Named numeric vector containing the parameter
#'               estimates from the parametric start.}
#'   \item{`logLik`}{The log-likelihood of the parametric starts. Is `NA`
#'               for the uniform start.}
#'
#' }
#'
#'
#' @details If `normalized` is `FALSE` and `start != "uniform"`, the resulting
#' density will not integrate to 1 in general.
#'
#' @seealso The `stats` package function [stats::density()].
#' @references
#'   Hjort, Nils Lid, and Ingrid K. Glad. "Nonparametric density estimation with a parametric start." The Annals of Statistics (1995): 882-904.
#'
#'   Jones, M. C., and D. A. Henderson. "Miscellanea kernel-type density estimation on the unit interval." Biometrika 94.4 (2007): 977-984.
#'
#'   Chen, Song Xi. "Probability density function estimation using gamma kernels." Annals of the Institute of Statistical Mathematics 52.3 (2000): 471-480.
#'
#'   Silverman, Bernard W. Density estimation for statistics and data analysis. Vol. 26. CRC press, 1986.
#'
#' @examples
#' ## Use gamma kernels to model positive data, the concentration of
#' ## theophylline
#'
#' concentration <- Theoph$conc + 0.001
#' plot(kdensity(concentration, start = "gamma", kernel = "gamma", adjust = 1 / 3),
#'   ylim = c(0, 0.15), lwd = 2, main = "Concentration of theophylline"
#' )
#' lines(kdensity(concentration, start = "gamma", kernel = "gaussian"),
#'   lty = 2, col = "grey", lwd = 2
#' )
#' lines(kdensity(concentration, start = "gaussian", kernel = "gaussian"),
#'   lty = 3, col = "blue", lwd = 2
#' )
#' lines(kdensity(concentration, start = "gaussian", kernel = "gamma", adjust = 1 / 3),
#'   lty = 4, col = "red", lwd = 2
#' )
#' rug(concentration)
#'
#' ## Using a density and and estimator from another package.
#'
#' skew_hyperbolic <- list(
#'   density   = SkewHyperbolic::dskewhyp,
#'   estimator = function(x) SkewHyperbolic::skewhypFit(x, printOut = FALSE)$param,
#'   support   = c(-Inf, Inf)
#' )
#'
#' kde <- kdensity(diff(LakeHuron), start = skew_hyperbolic)
#' plot(kde,
#'   lwd = 2, col = "blue",
#'   main = "Annual differences in water level (ft) of Lake Huron, 1875 - 1972"
#' )
#' lines(kde, plot_start = TRUE, lty = 2, lwd = 2) # Plots the skew hyperbolic density.
#' rug(diff(LakeHuron))
#'
#' kde$estimates # Also: coef(kde)
#' # Displays the parameter estimates:
#' #        mu     delta      beta        nu
#' # -1.140713  3.301112  2.551657 26.462469
#'
kdensity <- function(x, bw = NULL, adjust = 1, kernel = NULL, start = NULL,
                     support = NULL, na.rm = FALSE, normalized = TRUE,
                     tolerance = 0.01) {
  ## These tests are based purely on the data, and should be run no matter what.
  data.name <- deparse(substitute(x))
  has.na <- anyNA(x)

  assertthat::assert_that(!(has.na & !na.rm),
    msg = "x contains NAs and na.rm = FALSE."
  )

  x <- x[!is.na(x)] # This line is reached only if (has.na & !na.rm) is FALSE.

  ## 'kernel', 'start' and 'bw' can be custom made: In this case, they must
  ## be added to their environments.

  if (!is.null(start)) {
    if (!is.character(start)) {
      start_str <- deparse(substitute(start))
      add_start(start_str = start_str, start = start)
      start <- start_str
    }
  }

  if (!is.null(kernel)) {
    if (!is.character(kernel)) {
      kernel_str <- deparse(substitute(kernel))
      add_kernel(kernel_str = kernel_str, kernel = kernel)
      kernel <- kernel_str
    }
  }

  ## The case of bw == Inf is special! In this case, we return the parametric
  ## start itself.

  ## Now we handle the combinations of kernel, start and support.
  ## This is fancy defaults management.
  kss_list <- get_kernel_start_support(kernel, start, support)

  ## The strings are used for reporting and inside bandwidth functions,
  ## at must be preserved.
  start_str <- kss_list$start_str
  kernel_str <- kss_list$kernel_str

  ## We overwrite the kernel, start, and support with what we obtained from kss.
  kernel <- kss_list$kernel
  start <- kss_list$start
  support <- kss_list$support

  ## Tests for incompabibilities in the supplied values.
  support_compatible(kernel, start, support)

  ## We continue by checking if the supplied values make sense.
  if (!all(x <= support[2] & x >= support[1])) {
    stop(
      "The supplied data x is not contained in the support: (",
      support[1], ", ", support[2], ")."
    )
  }

  estimates <- start$estimator(x)
  parametric_start <- start$density

  # Name of the variable where the density is evaluated. Typically x.
  x_name <- names(formals(start$density))[1]

  parametric_start_vector <- function(x) {
    arguments <- list()
    arguments[[1]] <- x
    names(arguments)[1] <- x_name
    arguments <- append(arguments, as.list(estimates))
    do.call(parametric_start, arguments)
  }

  ## Takes care of the bandwidth. Can be either a double, a string, or a
  ## function taking the arguments data, kernel, start support.

  if (is.null(bw)) {
    bw <- get_standard_bw(kernel_str, start_str, support)
  }

  if (!is.numeric(bw)) {
    if (is.character(bw)) {
      bw_str <- bw
      bw <- get_bw(bw)(x, kernel_str, start_str, support)
    } else {
      bw_str <- deparse(substitute(bw))
      bw <- bw(x, kernel_str, start_str, support)
    }
  } else {
    bw_str <- "user supplied"
    if (bw == Inf) {
      msg <- "bw = Inf does not work with a uniform start."
      assertthat::assert_that(start_str != "uniform", start_str != "constant", msg = msg)
    }
  }

  ## The parameter h is computed. The basic bandwidth is h = bw*adjust for the
  ## normal kernel, and is adjusted for all the other kernels so that the sd
  ## of the kernel equals h.
  sd <- ifelse(!is.null(kss_list$kernel$sd), kss_list$kernel$sd, 1)
  h <- bw * adjust * sd

  ## The denominator can be computed.
  parametric_start_data <- parametric_start_vector(x)

  ## Now we handle normalization.

  normalization <- 1

  if (normalized & bw != Inf) {
    integrand <- function(y) {
      sapply(y, function(y) {
        mean(1 / h * kernel$kernel(y, x, h) /
          parametric_start_data) *
          parametric_start_vector(y)
      })
    }

    integral <- tryCatch(
      stats::integrate(integrand,
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

    assertthat::assert_that(!is.nan(integral$abs.error / integral$value),
      msg = msg
    )
    assertthat::assert_that(integral$abs.error / integral$value < tolerance,
      msg = msg
    )

    normalization <- integral$value
  }


  ## This is the main part of the returned object.

  return_function <- function(y) {
    if (missing(y)) {
      stop("A kdensity object is a density function. Call it with numerical input.")
    }

    if (h != Inf) {
      if (length(y) == 1) {
        mean(1 / h * kernel$kernel(y, x, h) * parametric_start_vector(y) /
          parametric_start_data) / normalization
      } else {
        sapply(y, function(y) {
          parametric_start_vector_y <- parametric_start_vector(y)
          1 / h * mean(kernel$kernel(y, x, h) * parametric_start_vector_y /
            parametric_start_data) / normalization
        })
      }
    } else {
      ## If bw = Inf, the parametric start is returned.
      arguments <- list()
      arguments[[1]] <- y
      names(arguments)[1] <- x_name
      arguments <- append(arguments, as.list(estimates))
      do.call(parametric_start, arguments)
    }
  }

  ## 'R6'-ish attributes:
  call <- match.call()
  range <- c(min(x), max(x))
  n <- length(x)
  logLik <- ifelse(start_str == "uniform" | start_str == "constant", NA,
    sum(parametric_start_vector(x))
  )

  ## S3 attributes:
  class(return_function) <- "kdensity"
  attr(return_function, "srcref") <- NULL
  return_function
}
