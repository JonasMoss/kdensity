### ===========================================================================
### GENERICS
### This package supports the following generics:
### -- plot, points, lines
### -- summary, print
### ===========================================================================

#' @export
`$.kdensity` = function(x, attr) {
  attr(x, attr)
}

#' @export
coef.kdensity = function(object, ...) object$estimate

#' @export
logLik.kdensity = function(object, ...) {
  msg = "'logLik' only makes sense for kdensity objects with a non-uniform parametric start."
  assertthat::assert_that(object$start != "uniform" & object$start != "constant", msg = msg)
  val = object$logLik
  attr(val, "nobs") = length(object$n)
  attr(val, "df")   = length(coef(object))
  class(val) = "logLik"
  val
}

#' @export
confint.kdensity = function(object, parm, level = 0.95, ...) {
  # Implement pointwise confidence intervals in some way.
}

#' Supplies a plotting range from a kdensity object.
#'
#' @param obj a kdensity object
#' @return a vector of size 1000, used for plotting.

get_range = function(obj) {
  support = attr(obj, "support")
  minimum = attr(obj, "range")[1]
  maximum = attr(obj, "range")[2]
  obj_range = maximum - minimum
  addition = obj_range/10
  xmin = max(minimum - addition, support[1])
  xmax = min(maximum + addition, support[2])
  seq(xmin, xmax, length.out = 1000)
}

#' Plot Method for Kernel Density Estimation
#'
#' The \code{plot} method for \code{kdensity} objects.
#'
#' @export
#' @param x a \code{kdensity} object.
#' @param range range of x values.
#' @param plot_start logical; if \code{TRUE}, plots the parametric start instead of the kernel density estimate.
#' @param zero_line logical; if \code{TRUE}, add a base line at \code{y = 0}.
#' @param ... further plotting parameters.
#' @return None.
#' @seealso \code{\link{kdensity}}
#' @examples
#' ## Using the data set "precip" to eye-ball the similarity between
#' ## a kernel fit, a parametric fit, and a kernel with parametric start fit.
#' kde_gamma = kdensity(precip, kernel = "gaussian", start = "gamma")
#' kde = kdensity(precip, kernel = "gaussian", start = "uniform")
#'
#' plot(kde_gamma, main = "Annual Precipitation in US Cities")
#' lines(kde_gamma, plot_start = TRUE, lty = 2)
#' lines(kde, lty = 3)
#' rug(precip)


plot.kdensity = function(x, range = NULL, plot_start = FALSE, zero_line = TRUE, ...) {
  plot_helper(x, range, plot_start, zero_line, ptype = "plot", ...)
}

#' @export
lines.kdensity = function(x, range = NULL, plot_start = FALSE, zero_line = TRUE, ...) {
  plot_helper(x, range, plot_start, zero_line, ptype = "lines", type = "l", ...)
}

#' @export
points.kdensity = function(x, range = NULL, plot_start = FALSE, zero_line = TRUE, ...) {
  plot_helper(x, range, plot_start, zero_line, ptype = "points", type = "p", ...)
}

#' Helper function for the plot methods.
#'
#' A helper function for the plot methods that does most of the work under
#' the hood.
#'
#' @param x A \code{kdensity} object.
#' @param range An optional range vector; like \code{x} in \code{plot.default}.
#' @param plot_start Logical; if \code{TRUE}, plots the parametric start only.
#' @param zero_line Logical; if \code{TRUE}, adds a line at \code{y = 0}.
#' @param ptype The kind of plot to make
#' @param ... Passed to plot.default.
#' @return None.

plot_helper = function(x, range = NULL, plot_start = FALSE, zero_line = TRUE, ptype = c("plot", "lines", "points"), ...) {

  ptype = match.arg(ptype)

  if(is.null(range)) range = get_range(x)

  ## Potential arguments in ellipses are handled here. They are modelled
  ## after the structure of the 'density' function.
  supplied = list(...)

  bw_string = NULL
  if(!is.null(x$bw_str)) bw_string = paste0(" ('", x$bw_str, "')")

  defaults = list(type = "l",
                  main = deparse(attr(x, "call")),
                  ylab = "Density",
                  xlab = paste0("N = ", x$n, "   Bandwidth = ", formatC(x$bw),
                               bw_string),
                  lwd  = 1)

  args = listmerge(x = defaults,
                   y = supplied)
  args$x = range

  if(plot_start) {
    start = attr(x, "start")

    msg = "To use 'plot_start = TRUE', supply a parametric start that is a proper density."
    assertthat::assert_that(!is.null(start), start != "uniform", msg = msg)

    start = get_start(start)

    parameters = attr(x, "estimates")
    parametric_start = start$density

    args$y = sapply(range, function(y) {
      do.call(parametric_start, as.list(c("x" = y, parameters)))})

  } else {
    args$y = x(range)
  }

  switch(ptype,
         plot   = do.call(graphics::plot, args),
         lines  = do.call(graphics::lines, args),
         points = do.call(graphics::points, args))

  if(zero_line) graphics::abline(h = 0, lwd = 0.1, col = "gray")

  invisible(x)

}

#' @export
print.kdensity <- function(x, ...)
{
  digits = list(...)$digits
  cat("\nCall:\n", deparse(attr(x, "call")), "\n\n",
      "Data:      ", attr(x, "data.name"), " (",attr(x, "n"), " obs.)\n",
      "Bandwidth: ", formatC(attr(x, "bw"), digits = digits), " ('", attr(x, "bw_str"), "')\n",
      "Support:   (", attr(x, "support")[1], ", ", attr(x, "support")[2],   ")\n",
      "Kernel:    ", attr(x, "kernel"), "\n",
      "Start:     ", attr(x, "start"), "\n\n",
      sep = "")
}

#' @export
summary.kdensity <- function(object, ...)
{
  digits = list(...)$digits
  parameters = attr(object, "estimates")
  params = NULL
  if(length(parameters) > 0)
  params = c("Parameter estimates:", "\n",
  sapply(1:length(parameters), function(i) paste0(names(parameters)[i], ": ", formatC(parameters[i], digits), "\n")),
  "\n")
  cat("\nCall: \n", deparse(attr(object, "call")), "\n\n",
      "Data:        ", attr(object, "data.name"), " (",attr(object, "n"), " obs.)\n",
      "Bandwidth:   ", formatC(attr(object, "bw"), digits = digits), " ('", attr(object, "bw_str"), "')\n",
      "Support:     (", attr(object, "support")[1], ", ", attr(object, "support")[2],   ")\n",
      "Kernel:      ", attr(object, "kernel"), "\n",
      "Start:       ", attr(object, "start"), "\n",
      "Range:       (", formatC(attr(object, "range")[1], digits), ", ", formatC(attr(object, "range")[2], digits),   ")\n",
      "NAs in data: ", attr(object, "has.na"), "\n",
      "Adjustment:  ", attr(object, "adjust"), "\n\n",
      params,
      sep = "")
}
