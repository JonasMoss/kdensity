#' Parametric starts
#'
#' A parametric start is a density function with associated estimator which
#' is used as a starting point in \code{kdensity}. Several parametric starts
#' are implemented, all with maximum likelihood estimation. Custom made
#' parametric starts are possible, see the Structure section.
#'
#' @section Built-in starts:
#' @usage NULL
#' @format NULL
#' @section Structure:
#'   The parametric start contains three elements: The density function, an
#'   estimation function, and the support of the density. The parameters of
#'   the density function must partially match the parameters of the estimator
#'   function. The estimator function takes one argument, a numeric vector,
#'   which is passed from \code{kdensity}.
#'
#' @examples
#'   start_exponential = list(
#'     density = dexp,
#'     estimator = function(data) {
#'       c(rate = 1/mean(data))
#'     },
#'     support   = c(0, Inf)
#'   )
#'
#'   start_inverse_gaussian = list(
#'     density = statmod::dinvgauss,
#'     estimator = function(data) {
#'       c(mean       = mean(data),
#'         dispersion = mean(1/data - 1/mean(data)))
#'     },
#'     support   = c(0, Inf)
#'   )
#'
#' @seealso \code{\link{kdensity}}; \code{\link{kernels}}; \code{\link{bandwidths}}
#'
#' @name starts

#' @rdname starts
#' @include builtin_starts_custom_maximum_likelihood.R
#' @usage NULL
#' @format NULL
#' @section Built-in starts:
#'    \code{uniform, constant}: Selecting the uniform start makes \code{kdensity}
#'    act like an ordinary kernel density estimator. The default value for any
#'    choice of kernel or support.
start_uniform = list(
  density   = function(x, dummy = NULL) 1,
  estimator = function(data) c(dummy = 1),
  support   = c(-Inf, Inf)
)

#' @rdname starts
#' @usage NULL
#' @format NULL
#' @section Built-in starts:
#'    \code{gaussian, normal}: The normal distribution. A natural choice for
#'    densities on R.
start_normal = list(
  density = dnorm,
  estimator = function(data) {
    c(mean = mean(data),
      sd   = sd(data))
  },
  support   = c(-Inf, Inf)
)

#' @rdname starts
#' @usage NULL
#' @format NULL
#' @section Built-in starts:
#'    \code{laplace}: The Laplace distribution.
start_laplace = list(
  density = function(x, mu, b) {
    1/(2*b)*exp(-1/b*abs(x-mu))
  },

  estimator = function(data) {
    c(mu = median(data),
      b  = mean(abs(data - mu)))
  },

  support   = c(-Inf, Inf)
)

#' @rdname starts
#' @usage NULL
#' @format NULL
#' @section Built-in starts:
#'    \code{exponential, gamma, lognormal, inverse_gaussian, weibull}: Densities
#'    supported on c(0, Inf).
start_exponential = list(
  density = dexp,
  estimator = function(data) {
    c(rate = 1/mean(data))
  },
  support   = c(0, Inf)
)

start_lognormal = list(
  density = dlnorm,
  estimator = function(data) {
    c(meanlog = mean(log(data)),
      sdlog   = sd(log(data)))
  },
  support   = c(0, Inf)
)

start_inverse_gaussian = list(
  density = statmod::dinvgauss,
  estimator = function(data) {
    c(mean       = mean(data),
      dispersion = mean(1/data - 1/mean(data)))
  },
  support   = c(0, Inf)
)

start_gamma = list(
  density   = dgamma,
  estimator = mlgamma,
  support   = c(0, Inf)
)

start_weibull = list(
  density   = dweibull,
  estimator = mlweibull,
  support   = c(0, Inf)
)

#' @rdname starts
#' @usage NULL
#' @format NULL
#' @section Built-in starts:
#'    \code{beta}: The beta distribution, supported on c(0, 1).
start_beta = list(
  density   = dbeta,
  estimator = mlbeta,
  support   = c(0, 1)
)
