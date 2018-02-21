#' Get densities and estimators from strings.
#'
#' @param start_str A string specifying the density of interest.
#' @return A list of two functions.

get_start = function(start_str) {

  assertthat::assert_that(is.character(start_str))

  if(start_str == "start_inverse_gaussian") {
    if(!("statmod" %in% rownames(utils::installed.packages()))) {
      stop("The option 'inverse_gaussian' requires the package 'statmod' to work.")
    }
  }

  parametric_start = switch(start_str,
    uniform          = start_uniform,
    normal           = start_normal,
    gaussian         = start_normal,
    gamma            = start_gamma,
    exponential      = start_exponential,
    inverse_gaussian = start_inverse_gaussian,
    lognormal        = start_lognormal,
    beta             = start_beta,
    kumar            = start_kumar,
    kumaraswamy      = start_kumaraswamy,
    laplace          = start_laplace,
    weibull          = start_weibull,
    gumbel           = start_gumbel,
    constant         = start_uniform,
    pareto           = start_pareto
  )

  if(is.null(parametric_start)) {
    if(exists(start_str)) {
      parametric_start = get(start_str)
    } else {
      stop(paste0("The supplied parametric start (",start_str,") is not implemented."))
    }
  }

  parametric_start

}

