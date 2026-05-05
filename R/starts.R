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
#'    20 built-in starts from the [univariateML][univariateML::univariateML-package] package, see
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

get_univariate_ml_support <- function(meta) {
  support <- meta$support

  if (isS4(support)) {
    return(list(
      bounds = support@.Data[1, ],
      type = support@type
    ))
  }

  list(
    bounds = support$bounds,
    type = support$type
  )
}

density_namespace_available <- function(density) {
  if (!grepl("::", density, fixed = TRUE)) {
    return(TRUE)
  }

  package_name <- sub("::.*", "", density)
  requireNamespace(package_name, quietly = TRUE)
}

get_density_and_support <- function(fun) {
  if (utils::packageVersion("univariateML") >= "1.5") {
    meta <- "univariateML::univariateML_metadata"
    meta <- eval(parser(paste0(meta,"[[paste0(\"ml\", fun)]]")))
    density <- meta$density
    support <- get_univariate_ml_support(meta)$bounds
    return(list(density = eval(parser(density)), support = support))
  }

  # nocov start
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
  # nocov end
}

if (utils::packageVersion("univariateML") >= "1.5") {
  meta <- "univariateML::univariateML_metadata"
  densities <- names(Filter(function(x) {
    get_univariate_ml_support(x)$type == "R" &&
      density_namespace_available(x$density)
  }, eval(parser(meta))))
  densities <- unname(sapply(densities, function(x) substring(x, 3)))
} else {
  densities <- univariateML::univariateML_models
}

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

## ---------------------------------------------------------------------------
## Accessors.
## ---------------------------------------------------------------------------

#' Get densities and estimators from strings.
#'
#' @keywords internal
#' @param start_str A string specifying the density of interest.
#' @return A list of two functions.

get_start <- function(start_str) {
  assert_(is.character(start_str))

  parametric_start <- starts_environment[[start_str]]

  msg <- paste0("The supplied parametric start ('", start_str, "') is not implemented.")
  assert_(!is.null(parametric_start), msg = msg)

  parametric_start
}

#' Add a new parametric start to `starts_environment`.
#'
#' @keywords internal
#' @param start_str A string giving the name of the density.
#' @param start The parametric start function.
#' @return None.

add_start <- function(start_str, start) {
  assert_(is.character(start_str))
  assert_(all(start_str == make.names(start_str)),
    msg = "The name of the parametric start is not valid. Use a short, valid name. (E.g. kdensity(x, start = gaussian), where gaussian is a predefined start function.)"
  )

  list_msg <- paste0("The parametric start ('", start_str, "') must be a list.")
  assert_(is.list(start), msg = list_msg)

  ## Checks for the right elements in start.
  density_msg <- paste0("The parametric start ('", start_str, "') must contain a function named 'density'.")
  estimator_msg <- paste0("The parametric start ('", start_str, "') must contain a function named 'estimator'.")
  support_msg <- paste0("The parametric start ('", start_str, "') must contain a vector named 'support'.")

  assert_(!is.null(start$density), msg = density_msg)
  assert_(!is.null(start$estimator), msg = estimator_msg)
  assert_(!is.null(start$support), msg = support_msg)

  assign(start_str, start, envir = starts_environment)
}
