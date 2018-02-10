### ===========================================================================
### GENERICS
### This package shoud support the following generics:
### -- plot, points, lines: Should work through the same interface.
### -- summary, print: Should be similar to 'density'
###
### Functions here are exported; the only other function is 'kdensity' itself.
### ===========================================================================

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
#' @param zero_lines logical; if \code{TRUE}, add a base line at \code{y = 0}.
#' @param ... further plotting parameters.
#' @return None.
#' @seealso \code{\link{kdensity}}
#' @examples
#' ## Using the data set "precip" to eye-ball the similarity between:
#' ## A kernel fit, a parametric fit, and a kernel with parametric start fit.
#' kde_gamma = kdensity(precip, kernel = "gaussian", start = "gamma")
#' kde = kdensity(precip, kernel = "gaussian", start = "uniform")
#'
#' plot(kde_gamma, main = "Annual Precipitation in US Cities")
#' lines(kde_gamma, plot_start = TRUE, lty = 2)
#' lines(kde, lty = 3)
#' rug(precip)


plot.kdensity = function(x, range = NULL, plot_start = FALSE, zero_line = TRUE, ...) {
  plot_helper(x, range, plot_start, zero_line, type = "plot", ...)
}

#' @export
lines.kdensity = function(x, range = NULL, plot_start = FALSE, zero_line = TRUE, ...) {
  plot_helper(x, range, plot_start, zero_line, type = "lines", ...)
}

#' @export
points.kdensity = function(x, range = NULL, plot_start = FALSE, zero_line = TRUE, ...) {
  plot_helper(x, range, plot_start, zero_line, type = "points", ...)
}

#' Helper function for the plot methods
#'
#' @inheritParams plot.kdensity
#' @param type the kind of plot to make
#' @return None.

plot_helper = function(x, range = NULL, plot_start = FALSE, zero_line = TRUE, type = c("plot", "lines", "points"), ...) {

  type = match.arg(type)

  if(is.null(range)) range = get_range(x)

  ## Potential arguments in ellipses are handled here. They are modelled
  ## after the structure of the 'density' function.
  supplied = list(...)
  defaults = list(type = "l",
                  main = deparse(attr(x, "call")),
                  ylab = "Density",
                  xlab = paste("N =", attr(x, "n"), "  Bandwidth =", formatC(attr(x, "bw"))),
                  lwd  = 1)

  args = listmerge(x = defaults,
                   y = supplied)
  args$x = range

  if(plot_start) {
    kernel = attr(x, "kernel")
    start = attr(x, "start")
    support = attr(x, "support")

    msg = "To use 'plot_start = TRUE', supply a parametric start that is a proper density."
    assertthat::assert_that(!is.null(start), start != "uniform", msg = msg)

    kss_list = get_kernel_start_support(kernel, start, support)
    start = kss_list$start
    start_str = kss_list$start_str

    parameters = attr(x, "estimates")
    parametric_start = start$density

    args$y = sapply(range, function(y) {
      do.call(parametric_start, as.list(c("x" = y, parameters)))})

  } else {
    args$y = x(range)
  }

  switch(type,
         plot   = do.call(plot, args),
         lines  = do.call(lines, args),
         points = do.call(points, args))

  if(zero_line) abline(h = 0, lwd = 0.1, col = "gray")

}

#' @export
print.kdensity <- function(obj, digits = NULL, ...)
{
  cat("\nCall:\n", deparse(attr(obj, "call")), "\n\n",
      "Data:      ", attr(obj, "data.name"), " (",attr(obj, "n"), " obs.)\n",
      "Bandwidth: ", formatC(attr(obj, "bw"), digits = digits), "\n",
      "Support:   (", attr(obj, "support")[1], ", ", attr(obj, "support")[2],   ")\n",
      "Kernel:    ", attr(obj, "kernel"), "\n",
      "Start:     ", attr(obj, "start"), "\n\n",
      sep = "")
}

#' @export
summary.kdensity <- function(obj, digits = NULL, ...)
{
  parameters = attr(obj, "estimates")
  cat("\nCall: \n", deparse(attr(obj, "call")), "\n\n",
      "Data:        ", attr(obj, "data.name"), " (",attr(obj, "n"), " obs.)\n",
      "Bandwidth:   ", formatC(attr(obj, "bw"), digits = digits), " ('", attr(obj, "bw_str"), "')\n",
      "Support:     (", attr(obj, "support")[1], ", ", attr(obj, "support")[2],   ")\n",
      "Kernel:      ", attr(obj, "kernel"), "\n",
      "Start:       ", attr(obj, "start"), "\n",
      "Range:       (", formatC(attr(obj, "range")[1], digits), ", ", formatC(attr(obj, "range")[2], digits),   ")\n",
      "NAs in data: ", attr(obj, "has.na"), "\n",
      "Adjustment:  ", attr(obj, "adjust"), "\n\n",
      "Parameter estimates:", "\n",
      sapply(1:length(parameters), function(i) paste0(names(parameters)[i], ": ", formatC(parameters[i], digits), "\n")),
      "\n",
      sep = "")
}
