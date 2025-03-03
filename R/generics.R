### ===========================================================================
### GENERICS
### This package supports the following generics:
### -- plot, points, lines
### -- summary, print
### ===========================================================================

#' @export
`[[.kdensity` <- function(x, i) {
  allowed_arg <- c(
    "x", "bw_str", "bw", "adjust", "h", "kernel_str", "kernel",
    "start", "start_str", "support", "data.name", "n", "range",
    "has.na", "na.rm", "normalized", "call", "estimates", "logLik"
  )

  i <- match.arg(i, allowed_arg)
  environment(x)[[i]]
}

#' @export
`[[<-.kdensity` <- function(x, i, value) {
  allowed_arg <- c(
    "x", "bw", "adjust", "kernel", "start", "support", "na.rm",
    "normalized"
  )

  i <- match.arg(i, allowed_arg)
  environment(x)$obj_name <- "x"
  args <- list(object = x)
  args[[i]] <- value
  do.call(update.kdensity, args)
  x
}


#' @export
`$<-.kdensity` <- function(x, name, value) {
  x[[name]] <- value
  x
}


#' @export
`$.kdensity` <- function(x, name) {
  x[[name]]
}

#' @export
update.kdensity <- function(object, ...) {
  current <- list(
    x = object$x,
    bw = object$bw_str,
    adjust = object$adjust,
    kernel = object$kernel_str,
    start = object$start_str,
    support = object$support,
    na.rm = object$na.rm,
    normalized = object$normalized
  )

  passed <- list(...)

  ## Part of a hack to make $<- and [[<- work.
  if (!is.null(environment(object)$obj_name)) {
    obj_name <- environment(object)$obj_name
  } else {
    obj_name <- deparse(substitute(object))
  }

  arg_names <- lapply(match.call(expand.dots = TRUE)[-1], deparse)
  args <- listmerge(current, passed, type = "template")
  new_object <- do.call(kdensity, args)

  if ("x" %in% names(arg_names)) {
    data.name <- arg_names$x
  } else {
    data.name <- object$data.name
  }

  call <- call("kdensity",
    x = data.name, adjust = args$adjust,
    kernel = args$kernel, start = args$start, support = args$support
  )

  environment(new_object)$call <- call

  environment(new_object)$data.name <- data.name
  assign(obj_name, new_object, envir = parent.frame())
}

#' @export
coef.kdensity <- function(object, ...) object$estimates

#' @export
logLik.kdensity <- function(object, ...) {
  msg <- "'logLik' only makes sense for kdensity objects with a non-uniform parametric start."
  assertthat::assert_that(object$start_str != "uniform" & object$start_str != "constant", msg = msg)
  val <- object$logLik
  attr(val, "nobs") <- length(object$n)
  attr(val, "df") <- length(stats::coef(object))
  class(val) <- "logLik"
  val
}

#' Supplies a plotting range from a kdensity object.
#'
#' @param obj A kdensity object.
#' @keywords internal
#' @return S vector of size 1000, used for plotting.

get_range <- function(obj) {
  support <- obj$support
  minimum <- obj$range[1]
  maximum <- obj$range[2]
  obj_range <- maximum - minimum
  addition <- obj_range / 10
  xmin <- max(minimum - addition, support[1])
  xmax <- min(maximum + addition, support[2])
  seq(xmin, xmax, length.out = 1000)
}

#' Plot, Lines and Points Methods for Kernel Density Estimation
#'
#' The `plot` method for `kdensity` objects.
#'
#' @export
#' @param x a `kdensity` object.
#' @param range range of x values.
#' @param plot_start logical; if `TRUE`, plots the parametric start instead of the kernel density estimate.
#' @param zero_line logical; if `TRUE`, add a base line at `y = 0`.
#' @param ... further plotting parameters.
#' @return None.
#' @seealso [kdensity()]
#' @examples
#' ## Using the data set "precip" to eye-ball the similarity between
#' ## a kernel fit, a parametric fit, and a kernel with parametric start fit.
#' kde_gamma <- kdensity(precip, kernel = "gaussian", start = "gamma")
#' kde <- kdensity(precip, kernel = "gaussian", start = "uniform")
#'
#' plot(kde_gamma, main = "Annual Precipitation in US Cities")
#' lines(kde_gamma, plot_start = TRUE, lty = 2)
#' lines(kde, lty = 3)
#' rug(precip)
plot.kdensity <- function(x, range = NULL, plot_start = FALSE, zero_line = TRUE, ...) {
  plot_helper(x, range, plot_start, zero_line, ptype = "plot", ...)
}

#' @rdname plot.kdensity
#' @export
lines.kdensity <- function(x, range = NULL, plot_start = FALSE, zero_line = TRUE, ...) {
  plot_helper(x, range, plot_start, zero_line, ptype = "lines", type = "l", ...)
}

#' @rdname plot.kdensity
#' @export
points.kdensity <- function(x, range = NULL, plot_start = FALSE, zero_line = TRUE, ...) {
  plot_helper(x, range, plot_start, zero_line, ptype = "points", type = "p", ...)
}

#' Helper function for the plot methods.
#'
#' A helper function for the plot methods that does most of the work under
#' the hood.
#'
#' @param x A `kdensity` object.
#' @param range An optional range vector; like `x` in `plot.default`.
#' @param plot_start Logical; if `TRUE`, plots the parametric start only.
#' @param zero_line Logical; if `TRUE`, adds a line at `y = 0`.
#' @param ptype The kind of plot to make
#' @param ... Passed to plot.default.
#' @keywords internal
#' @return None.

plot_helper <- function(x, range = NULL, plot_start = FALSE, zero_line = TRUE, ptype = c("plot", "lines", "points"), ...) {
  ptype <- match.arg(ptype)

  if (is.null(range)) range <- get_range(x)

  ## Potential arguments in ellipses are handled here. They are modelled
  ## after the structure of the 'density' function.
  supplied <- list(...)

  bw_string <- NULL
  if (!is.null(x$bw_str)) bw_string <- paste0(" ('", x$bw_str, "')")

  defaults <- list(
    type = "l",
    main = deparse(x$call),
    ylab = "Density",
    xlab = paste0(
      "N = ", x$n, "   Bandwidth = ", formatC(x$bw),
      bw_string
    ),
    lwd = 1
  )

  args <- listmerge(
    x = defaults,
    y = supplied
  )
  args$x <- range

  if (plot_start) {
    start <- x$start_str

    msg <- "To use 'plot_start = TRUE', supply a parametric start that is a proper density."
    assertthat::assert_that(!is.null(start), start != "uniform", msg = msg)

    start <- get_start(start)

    parameters <- stats::coef(x)
    parametric_start <- start$density

    args$y <- sapply(range, function(y) {
      do.call(parametric_start, as.list(c("x" = y, parameters)))
    })
  } else {
    args$y <- x(range)
  }

  switch(ptype,
    plot   = do.call(graphics::plot, args),
    lines  = do.call(graphics::lines, args),
    points = do.call(graphics::points, args)
  )

  if (zero_line) graphics::abline(h = 0, lwd = 0.1, col = "gray")

  invisible(x)
}

#' @export
print.kdensity <- function(x, ...) {
  digits <- list(...)$digits
  cat("\nCall:\n", deparse(x$call), "\n\n",
    "Data:      ", x$data.name, " (", x$n, " obs.)\n",
    "Bandwidth: ", formatC(x$bw, digits = digits), " ('", x$bw_str, "')\n",
    "Support:   (", x$support[1], ", ", x$support[2], ")\n",
    "Kernel:    ", x$kernel_str, "\n",
    "Start:     ", x$start_str, "\n\n",
    sep = ""
  )
  invisible(x)
}

#' @export
summary.kdensity <- function(object, ...) {
  digits <- list(...)$digits
  parameters <- object$estimates
  params <- NULL
  if (length(parameters) > 0) {
    params <- c(
      "Parameter estimates:", "\n",
      sapply(1:length(parameters), function(i) paste0(names(parameters)[i], ": ", formatC(parameters[i], digits), "\n")),
      "\n"
    )
  }
  cat("\nCall: \n", deparse(object$call), "\n\n",
    "Data:        ", object$data.name, " (", object$n, " obs.)\n",
    "Bandwidth:   ", formatC(object$bw, digits = digits), " ('", object$bw_str, "')\n",
    "Support:     (", object$support[1], ", ", object$support[2], ")\n",
    "Kernel:      ", object$kernel_str, "\n",
    "Start:       ", object$start_str, "\n",
    "Range:       (", formatC(object$range[1], digits), ", ", formatC(object$range[2], digits), ")\n",
    "NAs in data: ", object$has.na, "\n",
    "Adjustment:  ", object$adjust, "\n\n",
    params,
    sep = ""
  )
  invisible(object)
}
