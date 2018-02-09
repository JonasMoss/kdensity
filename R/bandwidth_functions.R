### ===========================================================================
### BANDWIDTH FUNCTIONS
###
### This file handles bandwidths. First we handle the conversion from strings
### to functions. The functions will take care of error handling. For instance,
### we will generate error if the kernel-start-support tupple is not
### compatible with the chosen bandwidth function. If it looks suspect, a
### warning is issued.
### ===========================================================================


#' Get bandwidth functions from string.
#'
#' @param bw a string specifying the density of interest.
#' @return a bandwidth function.
get_bw = function(bw) {
  switch(bw,
         "nrd0" = function(data, kernel_str, start_str, support) stats::bw.nrd0(data),
         "nrd"  = function(data, kernel_str, start_str, support) stats::bw.nrd(data),
         "ucv"  = function(data, kernel_str, start_str, support) stats::bw.ucv(data),
         "bcv"  = function(data, kernel_str, start_str, support) stats::bw.bcv(data),
         "SJ"   = function(data, kernel_str, start_str, support) stats::bw.SJ(data),
         "JH"   = bw.JH,
         "RHE"  = bw.RHE,
         stop("The supplied 'bw' is no among the supported alternatives.")
         )
}


#' Get a bandwidth string when 'bw' is unspecified.
#'
#' @param kernel_str supplied kernel,
#' @param start_str supplied parametric start
#' @param support supplied support
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

## Custom kernels.

bw.JH = function(data, kernel_str, start_str, support) {
  if(kernel_str != "gcopula") {
    warning("The bandwidth selection method JH is made for the asymmetric kernel 'gcopula'.")
  }

  if(support[1] < 0 | support[2] > 1) {
    warning("The bandwidth selection method JH is made for densities on the unit interval.")
  }

  ## The data is transfomed through qnorm, with singularities removed.
  transformed_data = qnorm(data)
  transformed_data = transformed_data[transformed_data != Inf & transformed_data != -Inf]
  sigma = sd(transformed_data)
  mu = mean(transformed_data)
  n = length(transformed_data)
  min(sigma * (2 * mu^2 * sigma^2 + 3*(1 - sigma^2)^2)^(-1/5)*n^(-1/5), 0.5)
}


bw.RHE = function(data, kernel_str, start_str, support) {
  assertthat::assert_that("EQL" %in% rownames(installed.packages()), msg =
    "The bandwidth function 'bw.RHE' requires the package 'EQL' to work.")
  max_degree = 5  # The maximum degree of the Hermite polynomials.
  n <- length(data)
  mu = mean(data)
  sigma = sd(data)
  z = (data - mu) / sigma

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

