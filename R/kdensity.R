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
#'   \item{`logLik`}{The log-likelihood of the *fitted parametric start*
#'               (not of the kernel density estimate): the sum of
#'               `log start(x_i)` at the maximum-likelihood estimates,
#'               evaluated directly on the log scale when supported by the
#'               start density. Is `NA` for the `uniform` and `constant`
#'               starts.}
#'
#' }
#'
#'
#' @details If `normalized` is `FALSE` and `start != "uniform"`, the resulting
#' density will not integrate to 1 in general.
#'
#' The `logLik()`, `AIC()`, and `BIC()` methods for `kdensity` objects describe
#' the *fitted parametric start* alone, not the kernel density estimate: they
#' use the start's maximum-likelihood log-likelihood, so they can be used to
#' compare candidate parametric starts. They are undefined for the `uniform` and
#' `constant` starts (`logLik()` returns `NA` and errors), which have no
#' parametric model.
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
kdensity <- function(x, bw = NULL, adjust = 1, kernel = NULL, start = NULL,
                     support = NULL, na.rm = FALSE, normalized = TRUE,
                     tolerance = 0.01) {
  ## These tests are based purely on the data, and should be run no matter what.
  data.name <- deparse(substitute(x))
  has.na <- anyNA(x)

  assert_(!(has.na & !na.rm),
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
  parametric_arguments <- c(
    stats::setNames(list(NULL), x_name),
    as.list(estimates)
  )

  parametric_start_vector <- function(x) {
    parametric_arguments[[1]] <- x
    do.call(parametric_start, parametric_arguments)
  }

  parametric_log_start_vector <- if (isTRUE(start$density_has_log)) {
    function(x) {
      arguments <- parametric_arguments
      arguments[[1]] <- x
      arguments$log <- TRUE
      do.call(parametric_start, arguments)
    }
  } else {
    function(x) log(parametric_start_vector(x))
  }

  ## Takes care of the bandwidth. Can be either a double, a string, or a
  ## function taking the arguments data, kernel, start support.

  bw_label <- deparse(substitute(bw))
  bw_info <- resolve_bandwidth(
    bw = bw,
    bw_label = bw_label,
    x = x,
    kernel_str = kernel_str,
    start_str = start_str,
    support = support
  )
  bw <- bw_info$bw
  bw_str <- bw_info$bw_str

  ## The parameter h is computed. The basic bandwidth is h = bw*adjust for the
  ## normal kernel, and is adjusted for all the other kernels so that the sd
  ## of the kernel equals h.
  sd <- if (!is.null(kss_list$kernel$sd)) kss_list$kernel$sd else 1
  h <- bw * adjust * sd

  ## The denominator can be computed.
  parametric_start_data <- parametric_start_vector(x)

  ## Now we handle normalization.

  normalization <- 1
  unnormalized_density <- NULL

  if (bw != Inf) {
    kernel_evaluator <- make_kernel_evaluator(
      x = x,
      h = h,
      kernel = kernel,
      parametric_start_data = parametric_start_data
    )

    unnormalized_density <- function(y) {
      kernel_evaluator(y) * parametric_start_vector(y)
    }

    if (normalized) {
      normalization <- resolve_normalization(
        unnormalized_density = unnormalized_density,
        start = start,
        kernel = kernel,
        support = support,
        h = h,
        x = x,
        estimates = estimates,
        tolerance = tolerance
      )
    }
  }


  ## This is the main part of the returned object.

  return_function <- function(y) {
    if (missing(y)) {
      stop("A kdensity object is a density function. Call it with numerical input.")
    }

    if (h != Inf) {
      unnormalized_density(y) / normalization
    } else {
      ## If bw = Inf, the parametric start is returned.
      parametric_start_vector(y)
    }
  }

  ## 'R6'-ish attributes:
  call <- match.call()
  range <- c(min(x), max(x))
  n <- length(x)
  logLik <- if (start_str == "uniform" || start_str == "constant") {
    NA
  } else {
    sum(parametric_log_start_vector(x))
  }

  ## S3 attributes:
  class(return_function) <- "kdensity"
  attr(return_function, "srcref") <- NULL
  return_function
}
