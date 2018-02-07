### ===========================================================================
### KERNEL FUNCTIONS
###
### This file handles the possible kernel functions. They fall into two
### categories. Standard symmetric kernel functions and non-standard assymetric
### kernel functions, such as the gamma kernel, beta kernel, and Gaussian
### copula kernel. Currently only symmetric kernel functions are supported.
### ===========================================================================

#' Helper function that gets a kernel function for kdensity.
#'
#' @param kernel a string specifying which kernel to use.
#' @return a kernel function of the format k(u) with integral normalized
#' to 1.

get_kernel = function(kernel) {
  switch(kernel,
         gaussian     = dnorm,
         laplace      = function(u) 1/2*exp(-abs(u)),
         epanechnikov = function(u) 3/4*(1-u^2)*(abs(u) <= 1),
         rectangular  = function(u) dunif(u, min = -1, max = 1),
         triangular   = function(u) (1-abs(u))*(abs(u) <= 1),
         biweight     = function(u) 15/16*(1-u^2)^2*(abs(u) <= 1),
         triweight    = function(u) 35/32*(1-u^2)^3*(abs(u) <= 1),
         tricube      = function(u) 70/81*(1-abs(u)^3)^3*(abs(u) <= 1),
         cosine       = function(u) (1+cos(pi*u))/2*(abs(u) <= 1),
         optcosine    = function(u) pi/4*cos(pi/2*u)*(abs(u) <= 1),
         uniform      = function(u) dunif(u, min = -1, max = 1)
  )
}

#' Helper function that gets a the kernel standard deviation.
#'
#' @param kernel a string specifying which kernel to use.
#' @return a kernel function of the format k(u) with standard deviation
#' normalized to 1.
#' @details stats::density use this concept in order to make different
#' kernels directly comparable. Most of the bandwidth selectors are based
#' one the Gaussian kernel, and translation requires the kernels to be
#' comparable.

get_kernel_sd = function(kernel) {
  switch(kernel,
         gaussian     = 1,
         laplace      = 1/sqrt(2),
         epanechnikov = sqrt(5),
         rectangular  = sqrt(3),
         triangular   = sqrt(6),
         biweight     = sqrt(7),
         triweight    = 3,
         tricube      = 3^(5/2)/sqrt(35),
         cosine       = 1/sqrt(1/3 - 2/pi^2),
         optcosine    = 1/sqrt(1-8/pi^2),
         uniform      = sqrt(3)
  )
}
