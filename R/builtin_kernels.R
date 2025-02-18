kernel_environment <- new.env(hash = FALSE)

#' Kernel functions
#'
#' Kernel functions are an important part of `kdensity`. This document
#' lists the available built-in functions and the structure of them. Any kernel
#' in the list can be used in `kdensity` by using `kernel = "kernel"`
#' for the intended kernel.
#'
#' Be careful combining kernels with compact support with parametric starts,
#' as the normalizing integral typically fails to converge. Use `gaussian`
#' instead.
#'
#' @section Symmetric kernels:
#' @section Asymmetric kernels:
#' @usage NULL
#' @format NULL
#' @section Structure:
#'    A kernel is a list containing two mandatory elements and one optional
#'    element. The mandatory element '`kernel`' is the kernel function.
#'    It takes arguments `y, x, h`, where `x` is the data supplied
#'    to `kdensity` and `y` is the point of evaluation. `h` is
#'    the bandwidth. Internally, the kernel function is evaluated as
#'    `1/h*kernel(y, x, h)`. It should be vectorized in `x`, but
#'    vectorization in `y` is not needed.
#'
#'    The second mandatory element is `support`, stating the domain of
#'    definition for the kernel. This is used to distinguish kernels on the
#'    unit interval / positive half-line from kernels on R.
#'
#'    `sd` is used for symmetric kernels, and states the standard error
#'    of the kernel. This is used to make kernels comparable to the Gaussian
#'    kernel when calculating bandwidths.
#' @examples
#' gaussian <- list(
#'   kernel = function(y, x, h) stats::dnorm((y - x) / h),
#'   sd = 1,
#'   support = c(-Inf, Inf)
#' )
#'
#' gcopula <- list(
#'   kernel = function(y, x, h) {
#'     rho <- 1 - h^2
#'     inside <- rho^2 * (qnorm(y)^2 + qnorm(x)^2) - 2 * rho * qnorm(y) * qnorm(x)
#'     exp(-inside / (2 * (1 - rho^2)))
#'   },
#'   support = c(0, 1)
#' )
#'
#' @seealso [kdensity()]; [parametric_starts()];
#' [bandwidths()].
#'
#' @section Symmetric kernels:
#'    `gaussian, normal`: The Gaussian kernel. The default argument when
#'    `starts` is supported on R.
#'    `epanechnikov, rectangular (uniform), triangular, biweight,
#'    cosine, optcosine`: Standard symmetric kernels, also used in
#'    [stats::density()].
#'    `tricube, triweight`: Standard symmetric kernels. Not supported by
#'    [stats::density()].
#'    `laplace`: Uses the Laplace density, also known as the double
#'    exponential density.
#' @section Asymmetric kernels:
#'    `gamma, gamma_biased`: The gamma kernel of Chen (2000). For use on the positive
#'    half-line. `gamma` is the recommended biased-corrected kernel.
#'    `gcopula`: The Gaussian copula kernel of Jones & Henderson (2007). For use
#'    on the unit interval.
#'    `beta, beta_biased`: The beta kernel of Chen (1999). For use on the unit interval.
#'    `beta` is the recommended bias-corrected kernel.
#' @references Chen, Song Xi. "Probability density function estimation using gamma kernels." Annals of the Institute of Statistical Mathematics 52.3 (2000): 471-480.
#' Jones, M. C., and D. A. Henderson. "Kernel-type density estimation on the unit interval." Biometrika 94.4 (2007): 977-984.
#' Chen, Song Xi. "Beta kernel estimators for density functions." Computational Statistics & Data Analysis 31.2 (1999): 131-145.
#' @name kernels
NULL

kernel_environment$gaussian <- list(
  kernel  = function(y, x, h) stats::dnorm((y - x) / h),
  sd      = 1,
  support = c(-Inf, Inf)
)

kernel_environment$normal <- kernel_environment$gaussian

kernel_environment$epanechnikov <- list(
  kernel = function(y, x, h) {
    u <- (x - y) / h
    3 / 4 * (1 - u^2) * (abs(u) <= 1)
  },
  sd = sqrt(5),
  support = c(-Inf, Inf)
)

kernel_environment$rectangular <- list(
  kernel  = function(y, x, h) dunif((x - y) / h, min = -1, max = 1),
  sd      = sqrt(3),
  support = c(-Inf, Inf)
)

kernel_environment$uniform <- kernel_environment$rectangular

kernel_environment$triangular <- list(
  kernel = function(y, x, h) {
    u <- (x - y) / h
    (1 - abs(u)) * (abs(u) <= 1)
  },
  sd = sqrt(6),
  support = c(-Inf, Inf)
)

kernel_environment$biweight <- list(
  kernel = function(y, x, h) {
    u <- (x - y) / h
    15 / 16 * (1 - u^2)^2 * (abs(u) <= 1)
  },
  sd = sqrt(7),
  support = c(-Inf, Inf)
)

kernel_environment$cosine <- list(
  kernel = function(y, x, h) {
    u <- (x - y) / h
    (1 + cos(pi * u)) / 2 * (abs(u) <= 1)
  },
  sd = 1 / sqrt(1 / 3 - 2 / pi^2),
  support = c(-Inf, Inf)
)

kernel_environment$optcosine <- list(
  kernel = function(y, x, h) {
    u <- (x - y) / h
    pi / 4 * cos(pi / 2 * u) * (abs(u) <= 1)
  },
  sd = 1 / sqrt(1 - 8 / pi^2),
  support = c(-Inf, Inf)
)

kernel_environment$triweight <- list(
  kernel = function(y, x, h) {
    u <- (x - y) / h
    35 / 32 * (1 - u^2)^3 * (abs(u) <= 1)
  },
  sd = 3,
  support = c(-Inf, Inf)
)

kernel_environment$tricube <- list(
  kernel = function(y, x, h) {
    u <- (x - y) / h
    70 / 81 * (1 - abs(u)^3)^3 * (abs(u) <= 1)
  },
  sd = 3^(5 / 2) / sqrt(35),
  support = c(-Inf, Inf)
)


kernel_environment$laplace <- list(
  kernel = function(y, x, h) {
    u <- (x - y) / h
    1 / 2 * exp(-abs(u))
  },
  sd = 1 / sqrt(2),
  support = c(-Inf, Inf)
)

kernel_environment$gcopula <- list(
  kernel = function(y, x, h) {
    rho <- 1 - h^2
    inside <- rho^2 * (qnorm(y)^2 + qnorm(x)^2) - 2 * rho * qnorm(y) * qnorm(x)
    exp(-inside / (2 * (1 - rho^2)))
  },
  support = c(0, 1)
)

kernel_environment$gamma <- list(
  kernel = function(y, x, h) {
    indices <- (y >= 2 * h)
    y_new <- y / h * indices + (1 / 4 * (y / h)^2 + 1) * (!indices)
    h * dgamma(x, y_new + 1, scale = h)
  },
  support = c(0, Inf)
)

kernel_environment$gamma_biased <- list(
  kernel = function(y, x, h) h * dgamma(x, y / h + 1, scale = h),
  support = c(0, Inf)
)

kernel_environment$beta <- list(
  kernel = function(y, x, h) {
    rho <- function(y, h) {
      2 * h^2 + 2.5 - sqrt(4 * h^4 + 6 * h^2 + 2.25 - y^2 - y / h)
    }

    N <- length(y)
    cut_points <- c(2 * h, 1 - 2 * h)
    sectioning <- findInterval(y, cut_points)

    y0 <- y[sectioning == 0]
    y1 <- y[sectioning == 1]
    y2 <- y[sectioning == 2]

    par1 <- c(rho(y0, h), y1 / h, y2 / h)
    par2 <- c((1 - y0) / h, (1 - y1) / h, rho(1 - y2, h))

    h * dbeta(x, par1, par2)
  },
  support = c(0, 1)
)

kernel_environment$beta_biased <- list(
  kernel  = function(y, x, h) h * dbeta(x, y / h + 1, (1 - y) / h + 1),
  support = c(0, 1)
)
