#' Helper function that gets a kernel function for kdensity.
#'
#' @param str a string specifying which kernel to use.
#' @return a kernel function of the format k(u) with integral normalized
#' to 1.

get_kernel = function(str) {

  assertthat::assert_that(is.character(str))

  kernel = switch(str,
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
         beta         = kernel_beta,
         beta_biased  = kernel_beta_biased
         )

  if(is.null(kernel)) {
    if(exists(str)) {
      kernel = get(str)
    } else {
      stop(paste0("The supplied kernel (",str,") is not implemented."))
    }
  }

  kernel

}
