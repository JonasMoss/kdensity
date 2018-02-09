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
         "nrd0" = function(data, kernel, start, support) stats::bw.nrd0(data),
         "nrd"  = function(data, kernel, start, support) stats::bw.nrd(data),
         "ucv"  = function(data, kernel, start, support) stats::bw.ucv(data),
         "bcv"  = function(data, kernel, start, support) stats::bw.bcv(data),
         "SJ"   = function(data, kernel, start, support) stats::bw.SJ(data)
         "JH"   = function(data, kernel, start, support) bw.JH

  )
}

## Custom kernels.

bw.JH = function(data, kernel_str, start_str, support) {
  if(kernel != "gcopula") {
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
