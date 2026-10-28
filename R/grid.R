#' Evaluate a kdensity fit on a fast grid
#'
#' `kdensity_grid()` evaluates a fitted [kdensity()] object on a grid. The
#' `"fft"` method uses linear binning and an FFT convolution for the built-in
#' Gaussian, rectangular, triangular, Epanechnikov, biweight, cosine, and
#' optimal-cosine kernels. It also supports the Gaussian-copula kernel through
#' an exact Gaussian-convolution representation on the probit scale.
#'
#' The FFT method approximates only the binning and interpolation steps. The
#' fitted parametric start and the normalization constant from `object` are
#' retained. Calling the `kdensity` object itself always uses exact evaluation.
#'
#' @param object A fitted `kdensity` object.
#' @param from,to Finite endpoints of the evaluation interval. By default,
#'   [get_range()] supplies a plotting interval.
#' @param n Number of returned grid points.
#' @param method One of `"auto"`, `"exact"`, or `"fft"`. `"auto"` uses the
#'   FFT backend for Gaussian and Gaussian-copula kernels and otherwise uses
#'   the exact evaluator. Compact kernels retain their faster exact sorted
#'   evaluator. `"fft"` explicitly requests FFT evaluation and reports an
#'   error when the fitted kernel is unsupported.
#' @param fft_n Size of the internal convolution grid. The default is at least
#'   four times `n` and at least 2048 points.
#' @param truncate Number of Gaussian standard deviations used to pad an FFT
#'   convolution. Compact kernels use their exact support instead.
#'
#' @return An object of class `"density"`. Its `method` element records whether
#'   exact or FFT evaluation was used.
#' @seealso [kdensity()]
#' @examples
#' set.seed(1)
#' fit <- kdensity(rbeta(1000, 2, 5), kernel = "gcopula", bw = 0.12)
#' fast <- kdensity_grid(fit, n = 512)
#' plot(fast)
#' @export
kdensity_grid <- function(object, from = NULL, to = NULL, n = 512L,
                          method = c("auto", "exact", "fft"),
                          fft_n = NULL, truncate = 8) {
  if (!inherits(object, "kdensity")) {
    stop("'object' must be a kdensity fit.")
  }
  method <- match.arg(method)
  n <- validate_grid_count(n, "n", 2L, "two")

  default_grid <- get_range(object)
  if (is.null(from)) {
    from <- min(default_grid)
  }
  if (is.null(to)) {
    to <- max(default_grid)
  }
  if (length(from) != 1L || length(to) != 1L ||
      !is.finite(from) || !is.finite(to) || from >= to) {
    stop("'from' and 'to' must be finite scalars with 'from' < 'to'.")
  }

  actual_method <- "exact"
  if (method != "exact") {
    support_problem <- fft_grid_support_problem(object, from, to)
    if (method == "auto") {
      kernel_id <- object$kernel$compiled_id
      actual_method <- if (
        is.null(support_problem) && kernel_id %in% c(1L, 11L)
      ) {
        "fft"
      } else {
        "exact"
      }
    } else if (!is.null(support_problem)) {
      stop("FFT grid evaluation is unavailable: ", support_problem)
    } else {
      actual_method <- "fft"
    }
  }

  if (actual_method == "fft") {
    if (!is.numeric(truncate) || length(truncate) != 1L ||
        is.na(truncate) || !is.finite(truncate) || truncate <= 0) {
      stop("'truncate' must be a single positive number.")
    }
    if (is.null(fft_n)) {
      fft_n <- max(2048, 4 * as.double(n))
    }
    fft_n <- validate_grid_count(fft_n, "fft_n", 512L, "512")
  } else {
    fft_n <- NA_integer_
  }

  grid <- make_kdensity_grid(
    object,
    from,
    to,
    n,
    truncate,
    use_probit = actual_method == "fft"
  )
  values <- if (actual_method == "fft") {
    evaluate_kdensity_fft_grid(
      object,
      grid = grid,
      fft_n = fft_n,
      truncate = truncate
    )
  } else {
    object(grid)
  }

  structure(
    list(
      x = grid,
      y = values,
      bw = object$bw,
      n = object$n,
      call = match.call(),
      data.name = object$data.name,
      has.na = object$has.na,
      method = actual_method,
      grid_n = n,
      fft_n = fft_n
    ),
    class = "density"
  )
}

validate_grid_count <- function(value, name, minimum, minimum_label) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < minimum || value != floor(value) ||
      value > .Machine$integer.max) {
    stop(
      "'", name, "' must be a single whole number of at least ",
      minimum_label, "."
    )
  }

  as.integer(value)
}

fft_grid_kernel_name <- function(kernel_id) {
  switch(as.character(kernel_id),
    "1" = "gaussian",
    "2" = "epanechnikov",
    "3" = "rectangular",
    "4" = "triangular",
    "5" = "biweight",
    "6" = "cosine",
    "7" = "optcosine",
    NULL
  )
}

fft_grid_support_problem <- function(object, from, to) {
  kernel_id <- object$kernel$compiled_id
  if (is.null(kernel_id) || !kernel_id %in% c(1:7, 11L)) {
    return("the fitted kernel has no FFT convolution representation.")
  }
  if (!is.finite(object$h) || object$h <= 0) {
    return("the fitted bandwidth is not finite and positive.")
  }

  environment_object <- environment(object)
  start_data <- environment_object$parametric_start_data
  normalization <- environment_object$normalization
  if (length(start_data) != object$n || any(!is.finite(start_data)) ||
      any(start_data <= 0)) {
    return("the fitted start has non-finite or non-positive data weights.")
  }
  if (length(normalization) != 1L || !is.finite(normalization) ||
      normalization <= 0) {
    return("the fitted normalization constant is not finite and positive.")
  }

  if (kernel_id == 11L) {
    if (from < 0 || to > 1) {
      return("Gaussian-copula FFT grids require [from, to] inside [0, 1].")
    }
    if (any(object$x <= 0 | object$x >= 1)) {
      return("Gaussian-copula FFT grids require data strictly inside (0, 1).")
    }
    rho <- 1 - object$h^2
    if (abs(rho) >= 1) {
      return("the Gaussian-copula bandwidth requires |1 - h^2| < 1.")
    }
  }

  NULL
}

make_kdensity_grid <- function(object, from, to, n, truncate,
                               use_probit = FALSE) {
  if (use_probit && identical(object$kernel$compiled_id, 11L) &&
      from >= 0 && to <= 1) {
    lower <- if (from == 0) -truncate else stats::qnorm(from)
    upper <- if (to == 1) truncate else stats::qnorm(to)
    grid <- stats::pnorm(seq(lower, upper, length.out = n))
    if (from == 0) {
      grid[[1L]] <- 0
    }
    if (to == 1) {
      grid[[n]] <- 1
    }
    return(grid)
  }
  seq(from, to, length.out = n)
}

# stats::density() gained `subdensity` and `ext` in R 4.4.0. On older R drop
# those two arguments and silence the deliberate sum(weights) != 1 warning; the
# grid values are unchanged (density self-extends its FFT grid regardless of
# `ext`), only this opt-in FFT path is affected.
density_fft <- function(...) {
  if (getRversion() >= "4.4.0") {
    return(stats::density(...))
  }
  # This compatibility path cannot run on the minimum R version used for
  # coverage, but remains necessary for older supported R releases.
  # nocov start
  args <- list(...)
  args$subdensity <- NULL
  args$ext <- NULL
  suppressWarnings(do.call(stats::density, args))
  # nocov end
}

evaluate_kdensity_fft_grid <- function(object, grid, fft_n, truncate) {
  environment_object <- environment(object)
  start_data <- environment_object$parametric_start_data
  start_density <- environment_object$parametric_start_vector
  normalization <- environment_object$normalization
  kernel_id <- object$kernel$compiled_id

  if (kernel_id == 11L) {
    rho <- 1 - object$h^2
    sigma <- sqrt(1 - rho^2)
    interior <- grid > 0 & grid < 1
    transformed <- rep(NA_real_, length(grid))
    transformed[interior] <- stats::qnorm(grid[interior])
    # The copula correction contains exp(z^2 / 2), which amplifies ordinary
    # FFT/interpolation roundoff in the extreme probit tails. At this limit the
    # amplification of machine epsilon reaches approximately sqrt(epsilon).
    # Exact evaluation of those few tail points prevents artificial boundary
    # spikes while retaining the fast convolution over the central grid.
    stable_limit <- sqrt(log(1 / .Machine$double.eps))
    fft_interior <- interior & abs(transformed) <= stable_limit
    if (sum(fft_interior) < 2L) {
      return(object(grid))
    }
    transformed_grid <- transformed[fft_interior]
    convolution <- density_fft(
      rho * stats::qnorm(object$x),
      weights = 1 / (object$n * object$h * start_data),
      bw = sigma,
      kernel = "gaussian",
      n = fft_n,
      from = min(transformed_grid),
      to = max(transformed_grid),
      cut = 0,
      ext = truncate,
      subdensity = TRUE
    )
    correction <- stats::approx(
      convolution$x,
      convolution$y,
      xout = transformed_grid,
      rule = 2
    )$y
    correction <- correction * sqrt(2 * pi) * object$h *
      exp(transformed_grid^2 / 2)
    values <- numeric(length(grid))
    values[fft_interior] <- correction * start_density(grid[fft_interior]) /
      normalization
    values[!fft_interior] <- object(grid[!fft_interior])
    return(values)
  } else {
    kernel_name <- fft_grid_kernel_name(kernel_id)
    kernel_sd <- object$kernel$sd
    convolution <- density_fft(
      object$x,
      weights = 1 / (object$n * start_data),
      bw = object$h / kernel_sd,
      kernel = kernel_name,
      n = fft_n,
      from = min(grid),
      to = max(grid),
      cut = 0,
      ext = if (kernel_id == 1L) {
        truncate
      } else {
        kernel_sd + sqrt(.Machine$double.eps)
      },
      subdensity = TRUE
    )
    correction <- stats::approx(
      convolution$x,
      convolution$y,
      xout = grid,
      rule = 2
    )$y
  }

  correction * start_density(grid) / normalization
}
