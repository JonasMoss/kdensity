#' Bandwidth Selectors
#'
#' Bandwidth selectors for \code{kdensity}. These are the available options
#' for the option \code{bw} in \code{kdensity}.
#'
#' @section Bandwidth selectors:
#'    \code{"nrd0", "nrd", "bcv", "SJ"}: Bandwidth selectors from \code{stats}.
#'    They are documented in \code{\link[stats]{bandwidth} stats:bandwidth}.
#'    "nrd0" is the standard bandwidth selector for symmetric kernels with
#'    constant parametric starts.
#'
#' @section Structure:
#'    The bandwidth selector is a function of four arguments: The data
#'    \code{x}, a kernel string \code{kernel}, a start string \code{start},
#'    and a support vector \code{support}. To obtain the functions associated
#'    with these strings, use \code{get_kernel} and \code{get_start}. The
#'    function should return a double.
#'
#' @seealso \code{\link{kdensity}}, \code{\link[stats]{bandwidth}} for the
#'    bandwidth selectors of \code{\link[stats]{density}}. In addition,
#'    \code{\link{kernels}}; \code{\link{starts}}
#'
#' @examples
#'    ## Not a serious bandwidth function.
#'    silly_width = function(x, kernel = NULL, start = NULL, support = NULL) {
#'      kernel = get_kernel(kernel)
#'      kernel$kernel(0.5)
#'    }
#' @name bandwidths

#' @rdname bandwidths
#' @usage NULL
#' @section Bandwidth selectors:
#'   \code{"JH"}: Selector for the Gaussian copula kernel, based on
#'   normal reference rule. Described in Jones & Henderson. The default method when
#'   the \code{gcopula} kernel is used in \code{kdensity}.
#' @references Jones, M. C., and D. A. Henderson. "Kernel-type density estimation on the unit interval." Biometrika 94.4 (2007): 977-984.
bw.JH = function(x, kernel = NULL, start = NULL, support = NULL) {

  # if(kernel != "gcopula") {
  #   warning("The bandwidth selection method JH is made for the asymmetric kernel 'gcopula'.")
  # }
  #
  # if(support[1] < 0 | support[2] > 1) {
  #   warning("The bandwidth selection method JH is made for densities on the unit interval.")
  # }

  ## The data is transfomed through qnorm, with singularities removed.
  transformed_x = stats::qnorm(x)
  transformed_x = transformed_x[transformed_x != Inf & transformed_x != -Inf]
  sigma = stats::sd(transformed_x)
  mu = mean(transformed_x)
  n = length(transformed_x)
  min(sigma * (2 * mu^2 * sigma^2 + 3*(1 - sigma^2)^2)^(-1/5)*n^(-1/5), 0.5)
}

#' @rdname bandwidths
#' @usage NULL
#' @section Bandwidth selectors:
#'   \code{"RHE"}: Selector for parametric starts with a symmetric kernel,
#'   based on a reference rule with Hermite polynomials.
#'   Described in Hjort & Glad (1995). The default method in \code{kdensity} when a parametric
#'   start is supplied and the kernel is symmetric.
#' @references Hjort, Nils Lid, and Ingrid K. Glad. "Nonparametric density estimation with a parametric start." The Annals of Statistics (1995): 882-904.
bw.RHE = function(x, kernel = NULL, start = NULL, support = NULL) {
  assertthat::assert_that("EQL" %in% rownames(utils::installed.packages()), msg =
                            "The bandwidth function 'bw.RHE' requires the package 'EQL' to work.")

  max_degree = 5  # The maximum degree of the Hermite polynomials.
  n <- length(x)
  mu = mean(x)
  sigma = stats::sd(x)
  z = (x - mu) / sigma

  ## Calculating the estimates of the robust Hermite polynomial coefficients.
  delta <- rep(0, max_degree)
  for (j in 2:max_degree) {
    hermite = EQL::hermite(sqrt(2) * z, j)
    delta[j] = mean(sqrt(2) * hermite * exp(-1/2 * z^2))
  }

  bw = (1/4)^(1/5) *
    (delta[2]^2 + delta[3]^2 + delta[4]^2/2 + delta[5]^2/6)^(-1/5) *
    sigma * n^(-1/5)
  return(bw)
}


#' @rdname bandwidths
#' @section Bandwidth selectors:
#'   \code{"ucv"}: Unbiased cross validation. The usual standard option for
#'   asymmetric kernels.

bw.ucv = function(x, kernel = NULL, start = NULL, support = NULL) {
  kernel_obj = get_kernel(kernel)
  start_obj = get_start(start)
  kernel_fun = kernel_obj$kernel
  start_density = start_obj$density
  start_estimator = start_obj$estimator
  start_support = start_obj$support
  n = length(x)

  # Name of the variable where the density is evaluated. Typically x.
  x_name = names(formals(start_density))[1]

  dstart = function(data, parameters) {
    sapply(data, function(datum) {
      arguments = as.list(c("x" = datum, parameters))
      names(arguments)[1] = x_name
      do.call(start_density, arguments)
    })
  }

  full_parameters = start_estimator(x)
  full_density_values = dstart(x, full_parameters)

  param_loo = list()
  for (i in 1:n) {
    param_loo[[i]] = start_estimator(x[-i])
  }
#
#   print(x)
#   print(kernel_fun)
#   print(start_density)
#   print(full_parameters)
#   print(support)

  parametric_start_vector = function(data) {
    sapply(data, function(datum) {
      arguments = as.list(c("x" = datum, full_parameters))
      names(arguments)[1] = x_name
      do.call(start_density, arguments)
    })
  }

  parametric_start_data = parametric_start_vector(x)

  obj_func = function(h) {
    kdensity_sqrd = kdensity_sq(x,
                                h = h,
                                kernel_fun = kernel_fun,
                                parametric_start_data = parametric_start_data,
                                parametric_start_vector = parametric_start_vector,
                                parametric_start = start_density,
                                support = support)

    term1 = stats::integrate(kdensity_sqrd,
                             lower = support[1],
                             upper = support[2])$value

    term2_vec = rep(0, n)
    for (i in 1:n) {
      term2_vec[i] = mean(kernel_fun(x, x[i], h) *
                          dstart(x, param_loo[[i]])) /
                          dstart(x[i], param_loo[[i]])
    }
    term2 = 2 * mean(term2_vec)

    obj_func_value = term1 - term2
    return(obj_func_value)
  }

  bw = optimize(obj_func, lower = 0.0001, upper = 10 * sd(x), tol = 0.0001)$minimum
  return(bw)
}


