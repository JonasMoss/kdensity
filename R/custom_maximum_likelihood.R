#' Fast maximum likelihood for the beta distribution.
#'
#' @param x The data, contained in the unit interval.
#' @param start Optional start point for the \code{nlm} function.
#' @param type Should we use a precomputed gradient function or Hessian?
#' @return a vector of two parameters, shape1 and shape2.
#' @details The option "none" is fastest at the moment. The starting value
#' is good.


mlbeta = function(x = NULL, start = NULL, type = c("none", "gradient", "full")) {
  type = match.arg(type, c("none","gradient","full"))

  if(!is.null(x)) {
    val1 = mean(log(x))
    val2 = mean(log(1-x))
  } else if (is.null(val1) & is.null(val2)) {
    stop("Either x or (val1 and val2) must be supplied.")
  }

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

  parameters = nlm(beta_objective, p = start, typsize = start)$estimate
  names(parameters) = c("shape1", "shape2")
  parameters

}
