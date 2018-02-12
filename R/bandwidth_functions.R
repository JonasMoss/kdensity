#' Bandwidth Selectors for Kernel Density Estimation
#'
#' Bandwidth selectors for \code{kdensity}.
#'
#' Bandwidth functions for parametric starts and asymmetric kernels currently include:
#'
#' \itemize{
#'   \item \strong{bw.JH:} Selector for the Gaussian copula kernel, based on
#'   normal reference rule. Described in Jones & Henderson. The default method when
#'   the \code{gcopula} kernel is used in \code{kdensity}.
#'   \item \strong{bw.RHE:} Selector for parametric starts with a symmetric kernel.
#'   Described in Hjort & Glad. The default method in \code{kdensity} when a parametric
#'   start is supplied and the kernel is symmetric.
#' }
#'
#' In addition to these,\code{kdensity} supports the bandwidth functions described in
#' \code{\link[stats]{bandwidth}}.
#'
#' @param x numeric vector containing the data.
#' @param kernel (optional) a \code{kernel} in list format.
#' @param start (optional) a \code{start} in list format.
#' @param support (optional) the support of the data.
#'
#' @seealso \code{\link{kdensity}}, \code{\link[stats]{bandwidth}} for the
#' bandwidth selectors for \code{\link[stats]{density}}.
#' @references
#' Hjort, Nils Lid, and Ingrid K. Glad. "Nonparametric density estimation with a parametric start." The Annals of Statistics (1995): 882-904.
#'
#' Jones, M. C., and D. A. Henderson. "Miscellanea kernel-type density estimation on the unit interval." Biometrika 94.4 (2007): 977-984.
#'
#' @name bandwidth_selectors

#' @rdname bandwidth_selectors
bw.JH = function(x, kernel = NULL, start = NULL, support = NULL) {
  if(kernel_str != "gcopula") {
    warning("The bandwidth selection method JH is made for the asymmetric kernel 'gcopula'.")
  }

  if(support[1] < 0 | support[2] > 1) {
    warning("The bandwidth selection method JH is made for densities on the unit interval.")
  }

  ## The data is transfomed through qnorm, with singularities removed.
  transformed_x = qnorm(x)
  transformed_x = transformed_x[transformed_x != Inf & transformed_x != -Inf]
  sigma = sd(transformed_x)
  mu = mean(transformed_x)
  n = length(transformed_x)
  min(sigma * (2 * mu^2 * sigma^2 + 3*(1 - sigma^2)^2)^(-1/5)*n^(-1/5), 0.5)
}

#' @rdname bandwidth_selectors
bw.RHE = function(x, kernel = NULL, start = NULL, support = NULL) {
  assertthat::assert_that("EQL" %in% rownames(installed.packages()), msg =
                            "The bandwidth function 'bw.RHE' requires the package 'EQL' to work.")
  max_degree = 5  # The maximum degree of the Hermite polynomials.
  n <- length(x)
  mu = mean(x)
  sigma = sd(x)
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

#' Get bandwidth functions from string.
#'
#' @param bw a string specifying the density of interest.
#' @return a bandwidth function.
get_bw = function(bw) {
  switch(bw,
         "nrd0" = function(data, kernel, start, support) stats::bw.nrd0(data),
         "nrd"  = function(data, kernel, start, support) stats::bw.nrd(data),
         "ucv"  = function(data, kernel, start, support) stats::bw.ucv(data),
         "bcv"  = function(data, kernel, start, support) stats::bw.bcv(data),
         "SJ"   = function(data, kernel, start, support) stats::bw.SJ(data),
         "JH"   = bw.JH,
         "RHE"  = bw.RHE,
         stop("The supplied 'bw' is no among the supported alternatives.")
         )
}


#' Get a bandwidth string when 'bw' is unspecified.
#'
#' @param kernel_str a kernel string
#' @param start_str a parametric start string.
#' @param support the support.
#' @return a bandwidth string.

get_standard_bw = function (kernel_str, start_str, support) {
  if(kernel_str == "gcopula") {
    bw = "JH"
  } else if (start_str != "uniform") {
    bw = "RHE"
  } else {
    bw ="nrd0"
  }
  bw
}
