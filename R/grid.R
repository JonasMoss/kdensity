#' Evaluate a kdensity fit on a fast grid
#'
#' `kdensity_grid()` evaluates a fitted [kdensity()] object on a grid. The
#' `"fft"` method uses weighted linear binning and an FFT convolution with the
#' fitted kernel, the algorithm of [stats::density()], for every built-in
#' symmetric kernel. It also supports the Gaussian-copula kernel through its
#' exact Gaussian-convolution representation on the probit scale.
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
#'   convolution. The Laplace kernel is padded by `truncate^2 / 2` bandwidths,
#'   where its tail is as small as the Gaussian tail at `truncate` standard
#'   deviations. Compact kernels use their exact support instead.
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

# Half-width of the window that contributes to an FFT grid value. Compact
# kernels (ids 2-9) vanish for |u| > 1, where u = (x - y) / h; the Gaussian is
# truncated at `truncate` standard deviations, and the Laplace kernel exp(-|u|)
# where it reaches the Gaussian tail level exp(-truncate^2 / 2).
fft_half_width <- function(kernel_id, h, truncate) {
  if (kernel_id == 1L) {
    return(truncate * h)
  }
  if (kernel_id == 10L) {
    return(truncate^2 / 2 * h)
  }
  h * (1 + sqrt(.Machine$double.eps))
}

fft_grid_support_problem <- function(object, from, to) {
  kernel_id <- object$kernel$compiled_id
  if (is.null(kernel_id) || !kernel_id %in% 1:11) {
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

# Weighted kernel sums  sum_i weights[i] * kernel_at(t - points[i])  at fft_n
# equally spaced nodes t covering [from - half_width, to + half_width], by
# linear binning (compiled, as in stats::density()) and a circular FFT
# convolution padded against wrap-around:
# the algorithm of stats::density(), with the caller's kernel. Points outside
# the node range lie more than half_width from [from, to] and are dropped.
# Roundoff can make sums slightly negative; they are clamped at zero.
binned_kernel_sum <- function(points, weights, kernel_at, half_width,
                              from, to, fft_n) {
  lower <- from - half_width
  upper <- to + half_width
  spacing <- (upper - lower) / (fft_n - 1)
  nodes <- lower + spacing * seq.int(0, fft_n - 1)

  binned <- linear_bin_cpp(points, weights, lower, spacing, fft_n)

  padded_n <- stats::nextn(2L * fft_n)
  lags <- spacing * seq.int(0, fft_n - 1)
  kernel_values <- numeric(padded_n)
  kernel_values[seq_len(fft_n)] <- kernel_at(lags)
  kernel_values[padded_n + 1L - seq_len(fft_n - 1L)] <- kernel_at(-lags[-1L])
  convolution <- stats::fft(
    stats::fft(c(binned, numeric(padded_n - fft_n))) *
      stats::fft(kernel_values),
    inverse = TRUE
  )
  list(
    x = nodes,
    y = pmax(Re(convolution[seq_len(fft_n)]) / padded_n, 0)
  )
}

evaluate_kdensity_fft_grid <- function(object, grid, fft_n, truncate) {
  environment_object <- environment(object)
  start_data <- environment_object$parametric_start_data
  start_density <- environment_object$parametric_start_vector
  normalization <- environment_object$normalization
  kernel_id <- object$kernel$compiled_id

  if (kernel_id == 11L) {
    h <- object$h
    rho <- 1 - h^2
    interior <- grid > 0 & grid < 1
    if (rho == 0 || sum(interior) < 2L) {
      # rho = 0 is the constant independence-copula kernel.
      return(object(grid))
    }
    # The copula density is symmetric in its arguments, so on the probit
    # scale z = qnorm(y) each term is a normal density in z:
    #   c(y, x_i) = dnorm(z, z_i / rho, scale) / (|rho| * dnorm(z_i)),
    # with z_i = qnorm(x_i) and scale = sqrt(1 - rho^2) / |rho|. The estimate
    # is a Gaussian convolution with no z-dependent factor, so the FFT's
    # roundoff is not amplified in the tails of the grid.
    transformed_grid <- stats::qnorm(grid[interior])
    transformed_x <- stats::qnorm(object$x)
    # Data deep in the probit tails carry weights 1 / dnorm(z_i) large enough
    # to swamp the roundoff of every other grid value; at this limit the
    # weights reach about 1 / sqrt(epsilon). Their few terms are summed
    # exactly instead.
    extreme <- abs(transformed_x) > sqrt(log(1 / .Machine$double.eps))
    scale <- sqrt(1 - rho^2) / abs(rho)
    regular <- !extreme
    convolution <- binned_kernel_sum(
      points = transformed_x[regular] / rho,
      weights = 1 / (object$n * start_data[regular] * abs(rho) *
        stats::dnorm(transformed_x[regular])),
      kernel_at = function(lag) stats::dnorm(lag, sd = scale),
      half_width = truncate * scale,
      from = min(transformed_grid),
      to = max(transformed_grid),
      fft_n = fft_n
    )
    sums <- stats::approx(
      convolution$x,
      convolution$y,
      xout = transformed_grid,
      rule = 2
    )$y + kernel_weighted_sum_cpp(
      y = grid[interior],
      x = object$x[extreme],
      weights = 1 / (object$n * h * start_data[extreme]),
      h = h,
      kernel_id = 11L,
      cache1 = transformed_x[extreme],
      cache2 = numeric()
    )
    values <- numeric(length(grid))
    values[interior] <- sums * start_density(grid[interior]) / normalization
    values[!interior] <- object(grid[!interior])
    return(values)
  } else {
    # Symmetric kernels depend on y - x only; evaluate the fitted kernel
    # itself at the lags so the FFT and exact paths share one definition.
    h <- object$h
    kernel_fun <- object$kernel$kernel
    convolution <- binned_kernel_sum(
      points = object$x,
      weights = 1 / (object$n * start_data),
      kernel_at = function(lag) kernel_fun(0, lag, h) / h,
      half_width = fft_half_width(kernel_id, h, truncate),
      from = min(grid),
      to = max(grid),
      fft_n = fft_n
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
