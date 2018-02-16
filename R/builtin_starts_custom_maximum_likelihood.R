#' Estimates the parameter of the Beta distribution using maximum likelihood
#'
#' Uses \code{stat::nlm} to estimate the parameters of the Beta distribution.
#'
#' @param x The data from which the estimate is to be computed.
#' @param start Optional starting parameter values for the minimization.
#' Passed to the \code{stats::nlm} function.
#' @param type Whether a dedicated \code{gradient} or \code{hessian} should be
#' passed to \code{stats::nlm}.
#' @return A named numeric vector with maximum likelihood estimates for
#' \code{shape1} and \code{shap2}.
#' @details For \code{type}, the option \code{none} is fastest. The algorithm
#' is a straight forwardimplementation of maximum likelihood for beta using
#' the sufficient statistics of the beta distribution,


mlbeta = function(x, start = NULL, type = c("none", "gradient", "hessian")) {
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

#' Estimates the parameter of the Gamma distribution using maximum likelihood
#'
#' Uses Newton-Raphson to estimate the parameters of the Gamma distribution.
#'
#' @param x The data from which the estimate is to be computed.
#' @param rel.tol Relative accuracy requested.
#' @param iterlim A positive integer specifying the maximum number of
#' iterations to be performed before the program is terminated.
#' @return A named numeric vector with maximum likelihood estimates for
#' \code{shape} and \code{rate}.
#' @references Choi, S. C, and R. Wette. "Maximum likelihood estimation of the parameters of the gamma distribution and their bias." Technometrics 11.4 (1969): 683-690.

mlgamma = function(x, rel_eps = 10^-10, iterlim = 100) {

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

  if(i == iterlim) {
    warning(paste0("The iteration limit (iterlim = ", iterlim, ") was reached",
                   " before the relative tolerance requirement (rel.tol = ",
                   rel.tol_str, ")."))
  }

  ## Given the shape, the rate is easy to compute.
  rate = shape/mean_hat

  c(shape = shape, rate = shape/mean_hat)

}


#' Estimates the parameter of a Weibull distribution by maximum likelihood
#'
#' Uses Newton-Raphson to estimate the parameters of the Weibull distribution.
#'
#' @param x The data from which the estimate is to be computed.
#' @param shape0 An optional starting value for the \code{shape} parameter.
#' @param rel.tol Relative accuracy requested.
#' @param iterlim A positive integer specifying the maximum number of
#' iterations to be performed before the program is terminated.
#'
#' @return A named numeric vector with maximum likelihood estimates for
#' \code{shape} and \code{scale}.
#' @export

mlweibull = function(x, initial = NULL, rel.tol = .Machine$double.eps^0.25,
                     iterlim = 100) {

  rel.tol_str = deparse(substitute(rel.tol))
  shape0 = 2
  log_x = log(x)
  l_hat = mean(log_x)
  log_xsq = log_x^2

  for(i in 1:iterlim) {
    shape0_lsum     = mean(x^shape0*log_x)
    shape0_lsum_sqr = mean(x^shape0*log_xsq)
    shape0_sum      = mean(x^shape0)
    A = shape0_lsum/shape0_sum
    B = shape0_lsum_sqr/shape0_sum
    top = 1/shape0 + l_hat - A
    bottom = -1/shape0^2 + A^2 - B
    shape = shape0 - top/bottom

    if(abs((shape0 - shape)/shape0) < rel.tol) break

    shape0 = shape
  }

  if(i == iterlim) {
    warning(paste0("The iteration limit (iterlim = ", iterlim, ") was reached",
                   " before the relative tolerance requirement (rel.tol = ",
                   rel.tol_str, ")."))
  }

  ## Given the shape, the scale is easy to compute.
  scale = (mean(x^shape))^(1/shape)
  c(shape = shape, scale = scale)
}
