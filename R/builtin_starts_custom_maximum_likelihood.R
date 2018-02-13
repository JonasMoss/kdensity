#' Fast maximum likelihood for the beta distribution.
#'
#' @param x The data, contained in the unit interval.
#' @param start Optional start point for the \code{nlm} function.
#' @param type Should we use a precomputed gradient function or Hessian?
#' @return a vector of two parameters, shape1 and shape2.
#' @details The option "none" is fastest at the moment. The algorithm is
#' a straight forward implementation of maximum likelihood for beta using
#' the sufficient statistics.


mlbeta = function(x, start = NULL, type = c("none", "gradient", "full")) {
  type = match.arg(type, c("none","gradient","full"))

  val1 = mean(log(x))
  val2 = mean(log(1-x))


  if(is.null(start)) {
    G1 = exp(val1)
    G2 = exp(val2)
    denom = 1/2*(1/(1-G1-G2))
    start = c(1/2 + G1*denom, 1/2 + G2*denom)
  }

  beta_objective = function(p) {
    lbeta(p[1],p[2]) - p[1]*val1 - p[2]*val2
  }

  if(type == "gradient" | type == "full") {
    attr(beta_objective, "gradient") = function(p) {
      digamma_alpha_beta = digamma(p[1] + p[2])
      c(digamma(p[1]) - digamma_alpha_beta - val1,
        digamma(p[2]) - digamma_alpha_beta - val2)
    }
  }

  if(type == "full") {
    attr(beta_objective, "hessian") = function(p) {
      trigamma_alpha_beta = -trigamma(p[1] + p[2])
      matrix(trigamma_alpha_beta, nrow = 2, ncol = 2) +
        diag(c(trigamma(p[1]), trigamma(p[2])))
    }}

  parameters = stats::nlm(beta_objective, p = start, typsize = start)$estimate
  names(parameters) = c("shape1", "shape2")
  parameters

}

#' Fast maximum likelihood for the gamma distribution
#'
#' @param x The data, contained in c(0, Inf)
#' @param rel_eps Relative epsilon comparison criterion.
#' @param max_iter Maximal number of Newton-Raphson iteratons.
#' @return A vector of two parameters, \code{shape} and \code{rate}.
#' @references Choi, S. C, and R. Wette. "Maximum likelihood estimation of the parameters of the gamma distribution and their bias." Technometrics 11.4 (1969): 683-690.


mlgamma = function(x, rel_eps = 10^-10, max_iter = 100) {

  mean_hat = mean(x)
  s = log(mean_hat) - mean(log(x))

  ## This start estimator is very close to the ML estimator already.
  shape = 1/(12*s)*(3 - s + sqrt((s-3)^2 + 24*s))

  ## The Newton-Raphson steps.
  for(i in 1:max_iter) {
    shape_next = shape - (1/(1/shape - trigamma(shape))*(log(shape) - digamma(shape) - s))
    if(abs((shape - shape_next)/shape) < rel_eps) break
    shape = shape_next
  }

  # gamma_objective = function(shape) {
  #   -((shape-1)*geom_hat - shape - shape*log_mean_hat + shape*log(shape) - lgamma(shape))
  # }
  #
  # optimize(gamma_objective, interval = c(0, 100))

  c(shape = shape, rate = shape/mean_hat)

}
