#' Parametric starts
#'
#' A parametric start is a density function with an associated estimator which
#' is used as a starting point in `kdensity`. Several parametric starts
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
#'   which is passed from `kdensity`.
#'
#' @examples start_exponential = list(
#'  density = stats::dexp,
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
#' @seealso [kdensity()]; [kernels()]; [bandwidths()]
#'
#' @section `kdensity` supports more than 20 built-in starts from the
#'    [univariateML]
#'    `uniform, constant`: Selecting the uniform start makes `kdensity`
#'    act like an ordinary kernel density estimator. The default value for any
#'    choice of kernel or support.
#'    `gaussian, normal`: The normal distribution. A natural choice for
#'    densities on the real line \eqn{(-\infty, \infty)}.
#'    `laplace, gumbel`: Distributions on  \eqn{(-\infty, \infty)}.
#'    `exponential, gamma, lognormal, inverse_gaussian, weibull`: Densities
#'    supported on the positive real line \eqn{(0, \infty)}.
#'    `beta, kumaraswamy`: The beta and Kumaraswamy distributions,
#'    supported on the unit interval \eqn{[0, 1]}.
#'    `pareto`: The Pareto distribution, supported on \eqn{[1, \infty)}.
#'    Has heavy tails.
#' @name parametric_starts
NULL

parser = function(str) parse(text = str)[[1]]

get_density_and_support = function(fun) {
  for(i in seq(length(body(fun)))) {

    if(length(body(fun)[[i]]) > 1) {
      if(body(fun)[[i]][[2]] == 'attr(object, "density")') {
        density = body(fun)[[i]][[3]]
      } else if (body(fun)[[i]][[2]] == 'attr(object, "support")') {
        support = body(fun)[[i]][[3]]
      }
    }
  }
  list(density = eval(parser(density)), support = support)
}

starts = new.env(hash = FALSE)

starts = lapply(univariateML::univariateML_models, function(name) {
  fun = eval(parser(paste0("univariateML::ml", name)))
  c(estimator = eval(parser(paste0("univariateML::ml",name))),
    get_density_and_support(fun))
})

names(starts) = univariateML::univariateML_models

## Some densities have variable supports, which is not supported yet.
starts$pareto = list(
  density   = function(x, alpha) alpha*x^(-alpha-1),
  estimator = function(x) 1/mean(log(x)),
  support   = c(1, Inf)
)

starts$power = NULL

## The uniform distribution is interpreted as uniform over the real line.
starts$unif= list(
  density   = function(x) rep(1, length(x)),
  estimator = function(data) NULL,
  support   = c(-Inf, Inf)
)

starts$constant = starts$unif
starts$uniform = starts$unif

## Aliases for densities.

starts$gaussian = starts$norm
starts$normal = starts$norm
starts$exponential = starts$exp
starts$lognormal = starts$lnorm
starts$inverse_gaussian = starts$invgauss
starts$wald = starts$invgauss

## Make starts_environments with evaled support.
for(i in seq_along(starts)) {
  starts[[i]]$support = eval(starts[[i]]$support)
}

starts_environment = as.environment(starts)

