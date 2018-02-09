### ===========================================================================
### SUPPORT HANDLERS
###
### The function in this file helps kdensity do its job: It assigns default
### supports to each parametric start. When non-symmetric kernels are
### implemented, these will affect the standard values as well.
### ===========================================================================

#' Helper function that generates automatic support.
#'
#' @param start the type of parametric start.
#' @param kernel the type of kernel.
#' @return a support proposal. Currently kernel does nothing.

get_support = function(start_str, kernel_str) {

  if(kernel_str == "gcopula") return(c(0, 1))

  switch(start_str,
         uniform          = c(-Inf, Inf),
         normal           = c(-Inf, Inf),
         gamma            = c(0, Inf),
         exponential      = c(0, Inf),
         inverse_gaussian = c(0, Inf),
         lognormal        = c(0, Inf),
         beta             = c(0, 1),
         laplace          = c(-Inf, Inf),
         c(-Inf, Inf)
  )
}
