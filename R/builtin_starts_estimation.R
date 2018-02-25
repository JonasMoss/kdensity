#' Estimates the parameter of the Beta distribution using maximum likelihood
#'
#' Uses \code{stat::nlm} to estimate the parameters of the Beta distribution.
#'
#' @keywords internal
#' @param x The data from which the estimate is to be computed.
#' @param start Optional starting parameter values for the minimization.
#' Passed to the \code{stats::nlm} function.
#' @param type Whether a dedicated \code{gradient} or \code{hessian} should be
#' passed to \code{stats::nlm}.
#' @return A named numeric vector with maximum likelihood estimates for
#' \code{shape1} and \code{shap2}.
#' @details For \code{type}, the option \code{none} is fastest.

mlbeta = function(x, start = NULL, type = c("none", "gradient", "hessian")) {
  type = match.arg(type)

  val1 = mean(log(x))
  val2 = mean(log(1-x))


  if(is.null(start)) {
    G1 = exp(val1)
    G2 = exp(val2)
    denom = 1/2*(1/(1-G1-G2))
    start = c(1/2 + G1*denom, 1/2 + G2*denom)
  }

  objective = function(p) {
    lbeta(p[1], p[2]) - p[1]*val1 - p[2]*val2
  }

  gradient = function(p) {
    digamma_alpha_beta = digamma(p[1] + p[2])
    c(digamma(p[1]) - digamma_alpha_beta - val1,
      digamma(p[2]) - digamma_alpha_beta - val2)
  }

  if(type == "gradient") {
    beta_objective = function(p) {
      result = objective(p)
      attr(result, "gradient") = gradient(p)
      result
    }
  } else if(type == "hessian") {
    hessian = function(p) {
      trigamma_alpha_beta = -trigamma(p[1] + p[2])
      matrix(trigamma_alpha_beta, nrow = 2, ncol = 2) +
        diag(c(trigamma(p[1]), trigamma(p[2])))
    }

    beta_objective = function(p) {
      result = objective(p)
      attr(result, "gradient") = gradient(p)
      attr(result, "hessian") = hessian(p)
      result
    }
  } else {
    beta_objective = objective
  }

  parameters = stats::nlm(beta_objective, p = start, typsize = start)$estimate
  names(parameters) = c("shape1", "shape2")
  parameters

}

#' Estimates the parameter of the Gamma distribution using maximum likelihood
#'
#' Uses Newton-Raphson to estimate the parameters of the Gamma distribution.
#'
#' @keywords internal
#' @param x The data from which the estimate is to be computed.
#' @param rel.tol Relative accuracy requested.
#' @param iterlim A positive integer specifying the maximum number of
#' iterations to be performed before the program is terminated.
#' @return A named numeric vector with maximum likelihood estimates for
#' \code{shape} and \code{rate}.
#' @references Choi, S. C, and R. Wette. "Maximum likelihood estimation of the parameters of the gamma distribution and their bias." Technometrics 11.4 (1969): 683-690.

mlgamma = function(x, rel.tol = .Machine$double.eps^0.25, iterlim = 100) {

  rel.tol_str = deparse(substitute(rel.tol))

  mean_hat = mean(x)
  s = log(mean_hat) - mean(log(x))

  ## This start estimator is very close to the ML estimator already.
  shape0 = 1/(12*s)*(3 - s + sqrt((s-3)^2 + 24*s))

  ## The Newton-Raphson steps.
  for(i in 1:iterlim) {
    shape = shape0 - (1/(1/shape0 - trigamma(shape0))*(log(shape0) - digamma(shape0) - s))
    if(abs((shape - shape0)/shape0) < rel.tol) break
    shape0 = shape
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
#' @keywords internal
#' @param x The data from which the estimate is to be computed.
#' @param shape0 An optional starting value for the \code{shape} parameter.
#' @param rel.tol Relative accuracy requested.
#' @param iterlim A positive integer specifying the maximum number of
#' iterations to be performed before the program is terminated.
#'
#' @return A named numeric vector with maximum likelihood estimates for
#' \code{shape} and \code{scale}.

mlweibull = function(x, shape0 = 2, rel.tol = .Machine$double.eps^0.25,
                     iterlim = 100) {

  rel.tol_str = deparse(substitute(rel.tol))
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

#' Estimates the parameter of a Gumbel distribution by maximum likelihood
#'
#' Uses Newton-Raphson to estimate the parameters of the Gumbel distribution.
#'
#' @keywords internal
#' @param x The data from which the estimate is to be computed.
#' @param scale0 An optional starting value for the \code{scale} parameter.
#' @param rel.tol Relative accuracy requested.
#' @param iterlim A positive integer specifying the maximum number of
#' iterations to be performed before the program is terminated.
#'
#' @return A named numeric vector with maximum likelihood estimates for
#' \code{shape} and \code{scale}.

mlgumbel = function(x, scale0 = 1, rel.tol = .Machine$double.eps^0.25,
                     iterlim = 100) {

  rel.tol_str = deparse(substitute(rel.tol))
  mean_x = mean(x)

  for(i in 1:iterlim) {

    A = sum(x*exp(-x/scale0))
    B = sum(exp(-x/scale0))
    C = sum(x^2*exp(-x/scale0))

    top = mean_x - scale0 - A/B
    bottom = -1 - 1/scale0^2*(C/B - (A/B)^2)

    scale = scale0 - top/bottom

    if(abs((scale0 - scale)/scale0) < rel.tol) break

    scale0 = scale
  }

  if(i == iterlim) {
    warning(paste0("The iteration limit (iterlim = ", iterlim, ") was reached",
                   " before the relative tolerance requirement (rel.tol = ",
                   rel.tol_str, ")."))
  }

  ## Given the shape, the scale is easy to compute.
  loc = -scale*log(mean(exp(-x/scale)))
  c(loc = loc, scale = scale)
}

#' Estimates the parameter of a Kumaraswamy distribution by maximum likelihood
#'
#' Uses Newton-Raphson to estimate the parameters of the Kumaraswamy distribution.
#'
#' @keywords internal
#' @param x The data from which the estimate is to be computed.
#' @param a0 An optional starting value for the \code{a} parameter.
#' @param rel.tol Relative accuracy requested.
#' @param iterlim A positive integer specifying the maximum number of
#'     iterations to be performed before the program is terminated.
#'
#' @return A named numeric vector with maximum likelihood estimates for
#'     \code{a} and \code{b}.
#'
#' @references Jones, M. C. "Kumaraswamy's distribution: A beta-type distribution with some tractability advantages." Statistical Methodology 6.1 (2009): 70-81.
#'
#'      Kumaraswamy, Ponnambalam. "A generalized probability density function for double-bounded random processes." Journal of Hydrology 46.1-2 (1980): 79-88.
#'

mlkumar = function(x, a0 = 1, rel.tol = .Machine$double.eps^0.25,
                    iterlim = 100) {

  rel.tol_str = deparse(substitute(rel.tol))

  logs = log(x)

  for(i in 1:iterlim) {

    xa = x^a0
    T1 = a0*mean(logs/(1-xa))
    T2 = a0*mean(logs*(xa/(1-xa)))
    T3 = mean(log(1-xa))
    f = 1/a0*(1 + T1 + T2/T3)

    C = mean(logs^2*(xa/(1-xa)^2))
    D = mean(logs^2*(xa/(1-xa)))

    T1diff = 1/a0*T1 + a0*C
    T2diff = 1/a0*T2  + 1/a0*T2^2 + a0*D

    fdiff = -1/a0^2*f + 1/a0*(T1diff + T2diff/T3 + 1/a0*(T2/T3)^2)

    a = a0 - f/fdiff
    if(abs((a0 - a)/a0) < rel.tol) break
    a0 = a

  }

  if(i == iterlim) {
    warning(paste0("The iteration limit (iterlim = ", iterlim, ") was reached",
                   " before the relative tolerance requirement (rel.tol = ",
                   rel.tol_str, ")."))
  }

  ## Given the shape, the scale is easy to compute.
  b = -1/mean(log(1 - x^a))
  c(a = a, b = b)
}
