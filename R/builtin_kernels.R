#' Kernel functions
#'
#' Kernel functions are an important part of \code{kdensity}. This document
#' lists the available built-in functions and the structure of them. Any kernel
#' in the list can be used in \code{kdensity} by using \code{kernel = "kernel"}
#' for the intended kernel.
#' @section Symmetric kernels:
#' @section Asymmetric kernels:
#' @usage NULL
#' @format NULL
#' @section Structure:
#'    A kernel is a list containing two mandatory elements and one optional
#'    element. The mandatory element '\code{kernel}' is the kernel function.
#'    It takes arguments \code{y, x, h}, where \code{x} is the data supplied
#'    to \code{kdensity} and \code{y} is the point of evaluation. \code{h} is
#'    the bandwidth. The kernel function is evaluated as
#'    \code{1/h*kerne(y, x, h)}, so multiply by \code{h} if that is needed. It
#'    should be vectorized in \code{x}, but vectorization in \code{y} is not
#'    needed.
#'
#'    The second mandatory element is \code{support}, stating the domain of
#'    definition for the kernel. This is used to distinguish kernels on the
#'    unit interval / positive half-line from kernels on R.
#'
#'    \code{sd} is used for symmetric kernels, and states the standard error
#'    of the kernel. This is used to make kernels comparable to the Gaussian
#'    kernel when calculating bandwidths.
#' @examples
#' gaussian = list(
#'   kernel  = function(y, x, h) dnorm((y-x)/h),
#'   sd = 1,
#'   support = c(-Inf, Inf)
#' )
#'
#' gcopula = list(
#'   kernel  = function(y, x, h) {
#'     rho = 1 - h^2
#'     inside = rho^2*(qnorm(y)^2 + qnorm(x)^2)-2*rho*qnorm(y)*qnorm(x)
#'     exp(-inside/(2*(1-rho^2)))
#'   },
#'   support = c(0, 1)
#' )
#'
#' @seealso \code{\link{kdensity}}; \code{\link{starts}}; \code{\link{bandwidths}}
#'
#' @name kernels

#' @rdname kernels
#' @usage NULL
#' @format NULL
#' @section Symmetric kernels:
#'    \code{gaussian}: The Gaussian kernel, defined as bla bla.
kernel_gaussian     = list(kernel  = function(y, x, h) dnorm((y-x)/h),
                           sd      = 1,
                           support = c(-Inf, Inf))

#' @rdname kernels
#' @usage NULL
#' @format NULL
#' @section Symmetric kernels:
#'    \code{laplace}: The Laplace kernel, defined as bla bla.
kernel_laplace      = list(kernel  = function(y, x, h) {
                                      u = (x - y)/h
                                      1/2*exp(-abs(u))
                                      },
                           sd      = 1/sqrt(2),
                           support = c(-Inf, Inf))

#' @rdname kernels
#' @usage NULL
#' @format NULL
#' @section Symmetric kernels:
#'    \code{epanechnikov}: The Epanechnikov kernel, defined as bla bla.
kernel_epanechnikov = list(kernel  = function(y, x, h) {
                                      u = (x - y)/h
                                      3/4*(1-u^2)*(abs(u) <= 1)
                                      },
                           sd      = sqrt(5),
                           support = c(-Inf, Inf))

#' @rdname kernels
#' @usage NULL
#' @format NULL
kernel_rectangular  = list(kernel  = function(y, x, h) dunif((x - y)/h, min = -1, max = 1),
                           sd      = sqrt(3),
                           support = c(-Inf, Inf))

#' @rdname kernels
#' @usage NULL
#' @format NULL
kernel_triangular   = list(kernel  = function(y, x, h) {
                                      u = (x - y)/h
                                      (1-abs(u))*(abs(u) <= 1)
                                     },
                           sd      = sqrt(6),
                           support = c(-Inf, Inf))

#' @rdname kernels
#' @usage NULL
#' @format NULL
kernel_biweight     = list(kernel  = function(y, x, h) {
                                      u = (x - y)/h
                                      15/16*(1-u^2)^2*(abs(u) <= 1)
                                     },
                           sd      = sqrt(7),
                           support = c(-Inf, Inf))

#' @rdname kernels
#' @usage NULL
#' @format NULL
kernel_triweight    = list(kernel  = function(y, x, h) {
                                      u = (x - y)/h
                                      35/32*(1-u^2)^3*(abs(u) <= 1)
                                     },
                           sd      = 3,
                           support = c(-Inf, Inf))

#' @rdname kernels
#' @usage NULL
#' @format NULL
kernel_tricube      = list(kernel = function(y, x, h) {
                                     u = (x - y)/h
                                     70/81*(1-abs(u)^3)^3*(abs(u) <= 1)
                                     },

                           sd      = 3^(5/2)/sqrt(35),

                           support = c(-Inf, Inf))

#' @rdname kernels
#' @usage NULL
#' @format NULL
kernel_cosine       = list(kernel  = function(y, x, h) {
                                      u = (x - y)/h
                                      (1+cos(pi*u))/2*(abs(u) <= 1)
                                     },
                           sd      = 1/sqrt(1/3 - 2/pi^2),
                           support = c(-Inf, Inf))

#' @rdname kernels
#' @usage NULL
#' @format NULL
kernel_optcosine    = list(kernel  = function(y, x, h) {
                                      u = (x - y)/h
                                      pi/4*cos(pi/2*u)*(abs(u) <= 1)
                                     },
                           sd      = 1/sqrt(1-8/pi^2),
                           support = c(-Inf, Inf)
                           )

#' @rdname kernels
#' @usage NULL
#' @format NULL
#' @section Asymmetric kernels:
#'    \code{gcopula}: The Gaussian copula kernel of Jones & Henderson (2007). For use
#'    on the unit interval.
#' @references Jones, M. C., and D. A. Henderson. "Kernel-type density estimation on the unit interval." Biometrika 94.4 (2007): 977-984.
kernel_gcopula = list(
  kernel= function(y, x, h) {
    rho = 1 - h^2
    inside = rho^2*(qnorm(y)^2 + qnorm(x)^2)-2*rho*qnorm(y)*qnorm(x)
    exp(-inside/(2*(1-rho^2)))
  },
  support = c(0, 1)
)

#' @rdname kernels
#' @usage NULL
#' @format NULL
#' @section Asymmetric kernels:
#'    \code{gamma, gamma_biased}: The gamma kernel of Chen (2000). For use on the positive
#'    half-line. \code{gamma} is the recommended biased-corrected kernel.
#' @references Chen, Song Xi. "Probability density function estimation using gamma kernels." Annals of the Institute of Statistical Mathematics 52.3 (2000): 471-480.
kernel_gamma = list(
  kernel  = function(y, x, h) {
    indices = (y >= 2*h)
    y_new = y/h*indices + (1/4*(y/h)^2 + 1)*(!indices)
    h*dgamma(x, y_new + 1, scale = h)
  },
  support = c(0, Inf)
)

#' @rdname kernels
#' @usage NULL
#' @format NULL
kernel_gamma_biased = list(
  kernel = function(y, x, h) h*dgamma(x, y/h + 1, scale = h),
  support = c(0, Inf)
)

#' @rdname kernels
#' @usage NULL
#' @format NULL
#' @section Asymmetric kernels:
#'    \code{beta, beta_biased}: The beta kernel of Chen (1999). For use on the unit interval.
#'    \code{beta} is the recommended bias-corrected kernel.
#' @references Chen, Song Xi. "Beta kernel estimators for density functions." Computational Statistics & Data Analysis 31.2 (1999): 131-145.
kernel_beta = list(
  kernel  = function(y, x, h) {

    rho = function(y, h) {
      2*h^2 + 2.5 - sqrt(4*h^4 + 6*h^2 + 2.25 - y^2 - y/h)
    }

    N = length(y)
    cut_points = c(2*h, 1 - 2*h)
    sectioning = findInterval(y, cut_points)

    y0 = y[sectioning == 0]
    y1 = y[sectioning == 1]
    y2 = y[sectioning == 2]

    par1 = c(rho(y0, h), y1/h, y2/h)
    par2 = c((1 - y0)/h, (1 - y1)/h, rho(1-y2, h))

    h*dbeta(x, par1, par2)
  }
  ,
  support = c(0, 1)
)

#' @rdname kernels
#' @usage NULL
#' @format NULL
kernel_beta_biased = list(
  kernel  = function(y, x, h) h*dbeta(x, y/h + 1, (1-y)/h + 1),
  support = c(0, 1)
)
