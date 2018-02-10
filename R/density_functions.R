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
#' @param start a string specifying the density of interest.
#' @param support a tuple containg the support. Used for error handling.
#' @return a list of two functions.

get_start = function(start) {
  if(start == "start_inverse_gaussian") {
    if(!("statmod" %in% rownames(installed.packages()))) {
      stop("The option 'inverse_gaussian' requires the package 'statmod' to work.")
    }
  }

  # ERROR HANDLING: Must be done elsewhere.
  # else if (start == "beta" & !all(support == c(0,1))) {
  #   stop("The option 'beta' requires support = c(0, 1).")
  # } else if (start %in% c("gamma, exponential, inverse_gaussian, lognormal") &!any(support < 0)) {
  #   stop("The supplied start option requires positive support.")
  # }

  switch(start,
         uniform          = start_uniform,
         normal           = start_normal,
         gamma            = start_gamma,
         exponential      = start_exponential,
         inverse_gaussian = start_inverse_gaussian,
         lognormal        = start_lognormal,
         beta             = start_beta,
         laplace          = start_laplace,
         start_uniform
         )
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
  support   = c(0, Inf)
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
  estimator = function(data) {
    shapes = nlm(function(par) {
      -mean(dgamma(data, par[1], par[2], log = TRUE))
    }, p = c(1,1))$estimate
    names(shapes) = c("shape", "rate")
    shapes
  },
  support   = c(0, Inf)
)
