#' Get bandwidth functions from string.
#'
#' @param bw a string specifying the density of interest.
#' @return a bandwidth function.
get_bw = function(bw) {

  assertthat::assert_that(is.character(bw))

  final_bw = switch(bw,
         "nrd0" = function(data, kernel, start, support) stats::bw.nrd0(data),
         "nrd"  = function(data, kernel, start, support) stats::bw.nrd(data),
         "bcv"  = function(data, kernel, start, support) stats::bw.bcv(data),
         "SJ"   = function(data, kernel, start, support) stats::bw.SJ(data),
         "ucv"  = bw.ucv,
         "JH"   = bw.JH,
         "RHE"  = bw.RHE
         )

  if(is.null(final_bw)) {
    if(exists(bw)) {
      final_bw = get(bw)
    } else {
      stop("The supplied 'bw' is no among the supported alternatives.")
    }
  }

  final_bw
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
