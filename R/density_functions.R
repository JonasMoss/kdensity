### ===========================================================================
### DENSITY FUNCTIONS
###
### This file handles the standard density functions. They are estimated
### efficiently by maximum likelihood. The rôle of kdensity_start_functions is
### to support handling of strings in the kdensity function, while the rest of
### file includes lists that follow the format of parametric start functions.
### ===========================================================================

#' Get densities and estimators from strings.
#'
#' @param start_str a string specifying the density of interest.
#' @return a list of two functions.

get_start = function(start_str) {

  assertthat::assert_that(is.character(start_str))

  if(start_str == "start_inverse_gaussian") {
    if(!("statmod" %in% rownames(installed.packages()))) {
      stop("The option 'inverse_gaussian' requires the package 'statmod' to work.")
    }
  }

  parametric_start = switch(start_str,
    uniform          = start_uniform,
    normal           = start_normal,
    gaussian         = start_normal,
    gamma            = start_gamma,
    exponential      = start_exponential,
    inverse_gaussian = start_inverse_gaussian,
    lognormal        = start_lognormal,
    beta             = start_beta,
    laplace          = start_laplace,
    constant         = start_uniform
  )

  if(is.null(parametric_start)) {
    if(exists(start_str)) {
      parametric_start = get(start_str)
    } else {
      stop(paste0("The supplied parametric start (",start_str,") is not implemented."))
    }
  }

  parametric_start

}


start_uniform = list(
  density   = function(x) 1,
  estimator = function(data) {},
  support   = c(-Inf, Inf)
)

start_normal = list(
  density = dnorm,
  estimator = function(data) {
    c(mean = mean(data),
      sd   = sd(data))
  },
  support   = c(-Inf, Inf)
)

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

start_beta = list(
  density   = dbeta,
  estimator = mlbeta,
  support   = c(0, 1)
)

start_gamma = list(
  density   = dgamma,
  estimator = mlgamma,
  support   = c(0, Inf)
)
