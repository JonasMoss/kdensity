starts_environment = new.env(hash = FALSE)

#' Parametric starts
#'
#' A parametric start is a density function with an associated estimator which
#' is used as a starting point in \code{kdensity}. Several parametric starts
#' are implemented, all with maximum likelihood estimation. Custom-made
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
#' @examples start_exponential = list(
#'  density = dexp,
#'  estimator = function(data) {
#'    c(rate = 1/mean(data))
#'  },
#'  support = c(0, Inf)
#' )
#'
#' start_inverse_gaussian = list(
#'  density = extraDistr::dwald,
#'  estimator = function(data) {
#'   c(mu = mean(data),
#'     lambda = mean(1/data - 1/mean(data)))
#'   },
#'  support = c(0, Inf)
#' )
#'
#' @seealso \code{\link{kdensity}}; \code{\link{kernels}}; \code{\link{bandwidths}}
#'
#' @section Built-in starts:
#'    \code{uniform, constant}: Selecting the uniform start makes \code{kdensity}
#'    act like an ordinary kernel density estimator. The default value for any
#'    choice of kernel or support.
#'    \code{gaussian, normal}: The normal distribution. A natural choice for
#'    densities on the real line \eqn{(-\infty, \infty)}.
#'    \code{laplace, gumbel}: Distributions on  \eqn{(-\infty, \infty)}.
#'    \code{exponential, gamma, lognormal, inverse_gaussian, weibull}: Densities
#'    supported on the positive real line \eqn{(0, \infty)}.
#'    \code{beta, kumaraswamy}: The beta and Kumaraswamy distributions,
#'    supported on the unit interval \eqn{[0, 1]}.
#'    \code{pareto}: The Pareto distribution, supported on \eqn{[1, \infty)}.
#'    Has heavy tails.
#' @name parametric_starts
NULL

starts_environment$uniform = list(
  density   = function(x) rep(1, length(x)),
  estimator = function(data) NULL,
  support   = c(-Inf, Inf)
)

starts_environment$constant = starts_environment$uniform

starts_environment$normal = list(
  density = dnorm,
  estimator = function(data) {
    c(mean = mean(data),
      sd   = sd(data))
  },
  support   = c(-Inf, Inf)
)

starts_environment$gaussian = starts_environment$normal

starts_environment$laplace = list(
  density = function(x, mu, b) {
    1/(2*b)*exp(-1/b*abs(x-mu))
  },

  estimator = function(data) {
    c(mu = median(data),
      b  = mean(abs(data - mu)))
  },

  support   = c(-Inf, Inf)
)

starts_environment$gumbel = list(
  density = function(x, loc, scale) {
      z = 1/scale*(x - loc)
      1/scale*exp(-(z + exp(-z)))
  },
  estimator = mlgumbel,
  support   = c(-Inf, Inf)
)

starts_environment$exponential = list(
  density = dexp,
  estimator = function(data) {
    c(rate = 1/mean(data))
  },
  support   = c(0, Inf)
)

starts_environment$lognormal = list(
  density = dlnorm,
  estimator = function(data) {
    c(meanlog = mean(log(data)),
      sdlog   = sd(log(data)))
  },
  support   = c(0, Inf)
)

if(requireNamespace("extraDistr", quietly = TRUE)) {
  starts_environment$inverse_gaussian = list(
    density = extraDistr::dwald,
    estimator = function(data) {
      c(mu       = mean(data),
        lambda   = mean(1/data - 1/mean(data)))
    },
    support   = c(0, Inf)
  )
} else {
  starts_environment$inverse_gaussian = list(
    density   = function(x) stop("Package 'extraDistr' required for option 'inverse_gaussian'."),
    estimator = function(x) stop("Package 'extraDistr' required for option 'inverse_gaussian'."),
    support   = c(0, Inf)
  )
}

starts_environment$wald = starts_environment$inverse_gaussian

starts_environment$gamma = list(
  density   = dgamma,
  estimator = mlgamma,
  support   = c(0, Inf)
)

starts_environment$weibull = list(
  density   = dweibull,
  estimator = mlweibull,
  support   = c(0, Inf)
  )

starts_environment$beta = list(
  density   = dbeta,
  estimator = mlbeta,
  support   = c(0, 1)
)

if(requireNamespace("extraDistr", quietly = TRUE)) {
  starts_environment$kumar = list(
    density   = extraDistr::dkumar,
    estimator = mlkumar,
    support   = c(0, 1)
  )
} else {
  starts_environment$kumar = list(
    density   = function(x) stop("Package 'extraDistr' required for option 'kumaraswsamy'."),
    estimator = function(x) stop("Package 'extraDistr' required for option 'kumaraswsamy'."),
    support   = c(0, Inf)
  )
}

starts_environment$pareto = list(
  density   = function(x, alpha) alpha*x^(-alpha-1),
  estimator = function(x) 1/mean(log(x)),
  support   = c(1, Inf)
)

