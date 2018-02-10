### ===========================================================================
### KERNEL FUNCTIONS
###
### This file handles the possible kernel functions. They fall into two
### categories. Standard symmetric kernel functions and non-standard assymetric
### kernel functions, such as the gamma kernel, beta kernel, and Gaussian
### copula kernel. Currently only symmetric kernel functions are supported.
###
### A user supplied kernel is a list containing a kernel function kernel that
### integrates to 1, and a float sd that equals the standard deviation of the
### kernel.
### ===========================================================================

#' Helper function that gets a kernel function for kdensity.
#'
#' @param kernel_str a string specifying which kernel to use.
#' @return a kernel function of the format k(u) with integral normalized
#' to 1.

get_kernel = function(kernel_str) {

  switch(kernel_str,
         gaussian     = kernel_gaussian,
         normal       = kernel_gaussian,
         laplace      = kernel_laplace,
         epanechnikov = kernel_epanechnikov,
         rectangular  = kernel_rectangular,
         triangular   = kernel_triangular,
         biweight     = kernel_biweight,
         triweight    = kernel_triweight,
         tricube      = kernel_tricube,
         cosine       = kernel_cosine,
         optcosine    = kernel_optcosine,
         uniform      = kernel_rectangular,
         gcopula      = kernel_gcopula,
         gamma        = kernel_gamma,
         gamma_biased = kernel_gamma_biased,
         stop(paste0("The supplied kernel (",kernel_str,") is not implemented."))
         )
}

## This is the list of pre-defined kernels.

kernel_gaussian     = list(kernel  = function(y, x, h) dnorm((y-x)/h),
                           sd      = 1,
                           support = c(-Inf, Inf))


kernel_laplace      = list(kernel  = function(y, x, h) {
                                      u = (x - y)/h
                                      1/2*exp(-abs(u))
                                      },
                           sd      = 1/sqrt(2),
                           support = c(-Inf, Inf))


kernel_epanechnikov = list(kernel  = function(y, x, h) {
                                      u = (x - y)/h
                                      3/4*(1-u^2)*(abs(u) <= 1)
                                      },
                           sd      = sqrt(5),
                           support = c(-Inf, Inf))


kernel_rectangular  = list(kernel  = function(y, x, h) dunif((x - y)/h, min = -1, max = 1),
                           sd      = sqrt(3),
                           support = c(-Inf, Inf))


kernel_triangular   = list(kernel  = function(y, x, h) {
                                      u = (x - y)/h
                                      (1-abs(u))*(abs(u) <= 1)
                                     },
                           sd      = sqrt(6),
                           support = c(-Inf, Inf))


kernel_biweight     = list(kernel  = function(y, x, h) {
                                      u = (x - y)/h
                                      15/16*(1-u^2)^2*(abs(u) <= 1)
                                     },
                           sd      = sqrt(7),
                           support = c(-Inf, Inf))


kernel_triweight    = list(kernel  = function(y, x, h) {
                                      u = (x - y)/h
                                      35/32*(1-u^2)^3*(abs(u) <= 1)
                                     },
                           sd      = 3,
                           support = c(-Inf, Inf))


kernel_tricube      = list(kernel = function(y, x, h) {
                                     u = (x - y)/h
                                     70/81*(1-abs(u)^3)^3*(abs(u) <= 1)
                                     },

                           sd      = 3^(5/2)/sqrt(35),

                           support = c(-Inf, Inf))

kernel_cosine       = list(kernel  = function(y, x, h) {
                                      u = (x - y)/h
                                      (1+cos(pi*u))/2*(abs(u) <= 1)
                                     },
                           sd      = 1/sqrt(1/3 - 2/pi^2),
                           support = c(-Inf, Inf))

kernel_optcosine    = list(kernel  = function(y, x, h) {
                                      u = (x - y)/h
                                      pi/4*cos(pi/2*u)*(abs(u) <= 1)
                                     },
                           sd      = 1/sqrt(1-8/pi^2),
                           support = c(-Inf, Inf)
                           )

kernel_gcopula      = list(kernel  = function(y, x, h) {
                                      rho = 1 - h^2
                                      inside = rho^2*(qnorm(y)^2 + qnorm(x)^2)-2*rho*qnorm(y)*qnorm(x)
                                      exp(-inside/(2*(1-rho^2)))
                                     },
                           sd      = 1,
                           support = c(0, 1))


kernel_gamma        = list(kernel  = function(y, x, h) {
                                      indices = (y >= 2*h)
                                      y_new = y/h*indices + (1/4*(y/h)^2 + 1)*(!indices)
                                      h*dgamma(x, y_new + 1, scale = h)
                                    },
                           sd      = 1,
                           support = c(0, Inf))

kernel_gamma_biased = list(kernel  = function(y, x, h) {
                                      h*dgamma(x, y/h + 1, scale = h)
                                    },
                           sd      = 1,
                           support = c(0, Inf))
