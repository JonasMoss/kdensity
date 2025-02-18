#' Parametric starts
#'
#' A parametric start is a density function with an associated estimator which
#' is used as a starting point in `kdensity`. Several parametric starts
#' are implemented, all with maximum likelihood estimation. Custom-made
#' parametric starts are possible, see the Structure section.
#'
#' @usage NULL
#' @format NULL
#' @section Structure:
#'   The parametric start contains three elements: The density function, an
#'   estimation function, and the support of the density. The parameters of
#'   the density function must partially match the parameters of the estimator
#'   function. The estimator function takes one argument, a numeric vector,
#'   which is passed from `kdensity`.
#'
#' @section Supported parametric starts: `kdensity` supports more than
#'    20 built-in starts from the [univariateML] package, see
#'    `univariateML::univariateML_models` for a list. Densities with variable
#'    support, `power`, are not supported. The `pareto` density has its
#'    support fixed to `(1,Inf)`. The
#'    options `uniform, constant` makes `kdensity` estimate a kernel
#'    density without parametric starts.
#' @examples start_exponential <- list(
#'   density = stats::dexp,
#'   estimator = function(data) {
#'     c(rate = 1 / mean(data))
#'   },
#'   support = c(0, Inf)
#' )
#'
#' start_inverse_gaussian <- list(
#'   density = extraDistr::dwald,
#'   estimator = function(data) {
#'     c(
#'       mu = mean(data),
#'       lambda = mean(1 / data - 1 / mean(data))
#'     )
#'   },
#'   support = c(0, Inf)
#' )
#'
#' @seealso [kdensity()]; [kernels()]; [bandwidths()]
#' @name parametric_starts
NULL

parser <- function(str) parse(text = str)[[1]]

get_density_and_support <- function(fun) {
  if (utils::packageVersion("univariateML") >= "1.5") {
    meta <- univariateML::univariateML_metadata[[paste0("ml", fun)]]
    density <- meta$density
    support <- meta$support@.Data[1, ]
    return(list(density = eval(parser(density)), support = support))
  }

  fun <- eval(parser(paste0("univariateML::ml", fun)))

  for (i in seq(length(body(fun)))) {
    if (length(body(fun)[[i]]) > 1) {
      if (body(fun)[[i]][[2]] == 'attr(object, "density")') {
        density <- body(fun)[[i]][[3]]
      } else if (body(fun)[[i]][[2]] == 'attr(object, "support")') {
        support <- body(fun)[[i]][[3]]
      }
    }
  }
  list(density = eval(parser(density)), support = support)
}

starts <- new.env(hash = FALSE)

densities <- names(Filter(\(x) x$support@type == "R", univariateML::univariateML_metadata))
densities <- unname(sapply(densities, \(x) substring(x, 3)))

starts <- lapply(densities, function(name) {
  c(
    estimator = eval(parser(paste0("univariateML::ml", name))),
    get_density_and_support(name)
  )
})

names(starts) <- densities

## Some densities have variable supports, which is not supported yet.
starts$pareto <- list(
  density   = function(x, alpha) alpha * x^(-alpha - 1),
  estimator = function(x) 1 / mean(log(x)),
  support   = c(1, Inf)
)

starts$power <- NULL

## The uniform distribution is interpreted as uniform over the real line.
starts$unif <- list(
  density   = function(x) rep(1, length(x)),
  estimator = function(data) NULL,
  support   = c(-Inf, Inf)
)

starts$constant <- starts$unif
starts$uniform <- starts$unif

## Aliases for densities.

starts$gaussian <- starts$norm
starts$normal <- starts$norm
starts$exponential <- starts$exp
starts$lognormal <- starts$lnorm
starts$inverse_gaussian <- starts$invgauss
starts$wald <- starts$invgauss

## Make starts_environments with evaled support.
for (i in seq_along(starts)) {
  starts[[i]]$support <- eval(starts[[i]]$support)
}

starts_environment <- as.environment(starts)
