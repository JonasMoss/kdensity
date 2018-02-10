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

#' Plotting generic for kdensity objects.
#' @export
#' @param obj a kdensity object
#' @param range the x-values the density is applied over (optional)
#' @param ... passed to plot.default.

plot.kdensity = function(obj, range = NULL, ...) {

  if(is.null(range)) range = get_range(obj)

  ## Potential arguments in ellipses are handled here. They are modelled
  ## after the structure of the 'density' function.
  supplied = list(...)
  defaults = list(type = "l",
                  main = deparse(attr(obj, "call")),
                  ylab = "Density",
                  xlab = paste("N =", attr(obj, "n"), "  Bandwidth =", formatC(attr(obj, "bw"))),
                  lwd  = 1)

  args = listmerge(x = defaults,
                   y = supplied)
  args$x = range
  args$y = obj(range)

  do.call(plot, args)

}

#' @export
lines.kdensity = function(obj, range = NULL, ...) {

  if(is.null(range)) range = get_range(obj)

  ## Potential arguments in ellipses are handled here. All go into the plot.
  supplied = list(...)
  defaults = list(type = "l",
                  lwd  = 1)

  args = listmerge(x = defaults,
                   y = supplied)
  args$x = range
  args$y = obj(range)

  do.call(lines, args)
}

#' @export
points.kdensity = function(obj, ...) {
  lines.kdensity(obj, type = "p", ...)
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
