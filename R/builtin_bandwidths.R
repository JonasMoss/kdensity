bw_environment <- new.env(hash = FALSE)

#' Bandwidth Selectors
#'
#' The available options for bandwidth selectors, passed as the `bw`
#' argument to `kdensity`.
#'
#' The bandwidth functions are not exported. They are members of the
#' environment `bw_environments`, and can be accessed by
#' `kdensity:::bw_environments`.
#'
#' @param x The input data.
#' @param kernel_str A string specifying the kernel, e.g. "gaussian."
#' @param start_str A string specifying the parametric start, e.g. "normal".
#' @param support The domain of definition for the kernel. (-Inf, Inf) for
#' symmetric kernels.
#'
#' @section Bandwidth selectors:
#'    `"nrd0", "nrd", "bcv", "SJ"`: Bandwidth selectors from `stats`.
#'    They are documented in `[bandwidth][stats::bandwidth] stats:bandwidth`.
#'    "nrd0" is the standard bandwidth selector for symmetric kernels with
#'    constant parametric starts.
#'
#'    `"ucv"`: Unbiased cross validation. The standard option for
#'    asymmetric kernels.
#'
#'    `"RHE"`: Selector for parametric starts with a symmetric kernel,
#'    based on a reference rule with Hermite polynomials.
#'    Described in Hjort & Glad (1995). The default method in `kdensity` when a parametric
#'    start is supplied and the kernel is symmetric.
#'
#'    `"JH"`: Selector for the Gaussian copula kernel, based on
#'    normal reference rule. Described in Jones & Henderson. The default method when
#'    the `gcopula` kernel is used in `kdensity`.
#'
#'
#' @section Structure:
#'    The bandwidth selector is a function of four arguments: The data
#'    `x`, a kernel string `kernel`, a start string `start`,
#'    and a support vector `support`. To obtain the functions associated
#'    with these strings, use `get_kernel` and `get_start`. The
#'    function should return a double.
#'
#' @seealso [kdensity()], [stats::bandwidth.kernel()] for the
#'    bandwidth selectors of [stats::density()]. In addition,
#'    [kernels()]; [parametric_starts()]
#'
#' @examples
#' ## Not a serious bandwidth function.
#' silly_width <- function(x, kernel = NULL, start = NULL, support = NULL) {
#'   rexp(1)
#' }
#' kdensity(mtcars$mpg, start = "gumbel", bw = silly_width)
#' @references
#' Jones, M. C., and D. A. Henderson. "Kernel-type density estimation on the unit interval." Biometrika 94.4 (2007): 977-984.
#' Hjort, Nils Lid, and Ingrid K. Glad. "Nonparametric density estimation with a parametric start." The Annals of Statistics (1995): 882-904.
#' @name bandwidths
NULL

bw_environment$JH <- function(x, kernel = NULL, start = NULL, support = NULL) {
  ## The data is transfomed through qnorm, with singularities removed.
  transformed_x <- stats::qnorm(x)
  transformed_x <- transformed_x[transformed_x != Inf & transformed_x != -Inf]
  sigma <- stats::sd(transformed_x)
  mu <- mean(transformed_x)
  n <- length(transformed_x)
  min(sigma * (2 * mu^2 * sigma^2 + 3 * (1 - sigma^2)^2)^(-1 / 5) * n^(-1 / 5), 0.5)
}

bw_environment$RHE <- function(x, kernel = NULL, start = NULL, support = NULL) {
  max_degree <- 5 # The maximum degree of the Hermite polynomials.
  n <- length(x)
  mu <- mean(x)
  sigma <- stats::sd(x)
  z <- (x - mu) / sigma

  ## Calculating the estimates of the robust Hermite polynomial coefficients.
  delta <- rep(0, max_degree)
  for (j in 2:max_degree) {
    hermite <- EQL::hermite(sqrt(2) * z, j)
    delta[j] <- mean(sqrt(2) * hermite * exp(-1 / 2 * z^2))
  }

  bw <- (1 / 4)^(1 / 5) *
    (delta[2]^2 + delta[3]^2 + delta[4]^2 / 2 + delta[5]^2 / 6)^(-1 / 5) *
    sigma * n^(-1 / 5)
  return(bw)
}

bw_environment$nrd0 <- function(data, kernel, start, support) stats::bw.nrd0(data)
bw_environment$nrd <- function(data, kernel, start, support) stats::bw.nrd(data)
bw_environment$bcv <- function(data, kernel, start, support) stats::bw.bcv(data)
bw_environment$SJ <- function(data, kernel, start, support) stats::bw.SJ(data)

bw_environment$ucv <- function(x, kernel = NULL, start = NULL, support = NULL) {
  ## We check for the combination start == "uniform" and kernel == "gaussian",
  ## as this is handled by stats::density's related functions.

  if (!is.null(start)) {
    if (start == "constant" | start == "uniform") {
      if (!is.null(get_kernel(kernel)$sd)) {
        return(stats::bw.ucv(x))
      }
    }
  }


  ## If we are here, we must do our own cross-validation.
  kernel_obj <- get_kernel(kernel)
  start_obj <- get_start(start)
  kernel_fun <- kernel_obj$kernel
  start_density <- start_obj$density
  start_estimator <- start_obj$estimator
  start_support <- start_obj$support
  n <- length(x)

  # Name of the variable where the density is evaluated. Typically x.
  x_name <- names(formals(start_density))[1]

  full_parameters <- start_estimator(x)

  arguments <- list()
  arguments[[1]] <- x
  names(arguments)[1] <- x_name
  arguments <- append(arguments, as.list(full_parameters))

  dstart <- function(data, parameters) {
    arguments[[1]] <- data
    if (length(parameters) > 0) {
      for (i in 1:length(parameters)) arguments[[i + 1]] <- parameters[[i]]
    }
    do.call(start_density, arguments)
  }

  param_loo <- lapply(1:n, function(i) start_estimator(x[-i]))

  parametric_start_vector <- function(data) dstart(data, full_parameters)

  parametric_start_data <- parametric_start_vector(x)

  obj_func <- function(h) {
    kdensity_sqrd <- kdensity_sq(x,
      h = h,
      kernel_fun = kernel_fun,
      parametric_start_data = parametric_start_data,
      parametric_start_vector = parametric_start_vector,
      parametric_start = start_density,
      support = support
    )

    term1 <- stats::integrate(kdensity_sqrd,
      lower = support[1],
      upper = support[2]
    )$value

    term2_vec <- rep(0, n)
    for (i in 1:n) {
      # Will not work fro asymmetric? Difference between x and y in kernel.
      term2_vec[i] <- mean(1 / h * kernel_fun(x[-i], x[i], h) /
        dstart(x[-i], param_loo[[i]])) *
        dstart(x[i], param_loo[[i]])
    }
    term2 <- 2 * mean(term2_vec)

    obj_func_value <- term1 - term2
    return(obj_func_value)
  }


  # The range of allowable bandwidths vary from kernel to kernel.

  eps <- 10^-10

  if (kernel == "gcopula" | kernel == "beta" | kernel == "beta_biased") {
    using_str <- "JH"
    using <- bw_environment$JH(x)
    lower <- 1 / 4 * using
    upper <- 1 / 4 - eps
  } else if (start == "constant" | start == "uniform") {
    using_str <- "nrd0"
    using <- stats::bw.nrd0(x)
    lower <- 1 / 5 * using
    upper <- 5 * using
  } else {
    using_str <- "RHE"
    using <- bw_environment$RHE(x)
    lower <- 1 / 5 * using
    upper <- 5 * using
  }

  bw <- tryCatch(stats::optimize(obj_func, lower = lower, upper = upper, tol = 0.0001)$minimum,
    error = function(e) {
      warning(paste0("Integration failed when finding bandwidth using 'ucv'. Using '", using_str, "' instead."), call. = FALSE)
      using
    }
  )
  return(bw)
}
