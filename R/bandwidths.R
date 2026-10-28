bw_environment <- new.env(hash = FALSE)

#' Bandwidth Selectors
#'
#' The available options for bandwidth selectors, passed as the `bw`
#' argument to `kdensity`.
#'
#' The bandwidth functions are not exported. They are members of the
#' environment `bw_environments`, and can be accessed by
#' `kdensity:::bw_environments`.
#'
#' @param x The input data.
#' @param kernel_str A string specifying the kernel, e.g. "gaussian."
#' @param start_str A string specifying the parametric start, e.g. "normal".
#' @param support The domain of definition for the kernel. (-Inf, Inf) for
#' symmetric kernels.
#'
#' @section Defaults:
#'    Every default is a closed-form rule, so no default requires
#'    cross-validation:
#'
#'    | Kernel | Uniform start | Parametric start |
#'    |---|---|---|
#'    | symmetric (`gaussian`, ...) | `"nrd0"` | `"RHE"` |
#'    | `gamma`, `gamma_biased` | `"GR"` | `"GR"` |
#'    | `beta`, `beta_biased` | `"HS"` | `"HS"` |
#'    | `gcopula` | `"SJ"` | `"SJ"` |
#'    | custom kernel without `sd` | `"ucv"` | `"ucv"` |
#'
#'    The rules use only the data, not the fitted parametric start. When the
#'    start fits well, the correction it leaves for the kernel is smoother than
#'    the density itself, so a larger bandwidth would be better; the rules then
#'    undersmooth slightly. `"ucv"` takes the start into account at the cost of
#'    a computation that grows with the square of the sample size.
#'
#' @section Bandwidth selectors:
#'    `"nrd0", "nrd", "bcv"`: Bandwidth selectors from `stats`, documented
#'    in [stats::bw.nrd0()]. They are on the scale of the Gaussian kernel's
#'    standard deviation; `kdensity` rescales them for the other symmetric
#'    kernels.
#'
#'    `"SJ"`: The Sheather & Jones (1991) plug-in selector [stats::bw.SJ()].
#'    For the `gcopula` kernel it is applied on the probit scale: there the
#'    Gaussian copula estimator is a Gaussian kernel estimate with bandwidth
#'    `s = sqrt(1 - rho^2)`, where `rho = 1 - h^2`, so the selected `s` is
#'    converted to `h = sqrt(1 - sqrt(1 - s^2))`. If the selector fails on
#'    degenerate data, `"JH"` is used with a warning.
#'
#'    `"ucv"`: Unbiased cross-validation. With a parametric start it refits the
#'    start without each observation, so its cost grows with the square of the
#'    sample size.
#'
#'    `"HS"`: Hallberg Szabadváry's closed-form reference rule for the
#'    bias-corrected beta kernel, based on a beta reference density fitted by
#'    moments. Used for both beta kernels, with a fallback heuristic when the
#'    rule's moment conditions fail.
#'
#'    `"GR"`: Chen's AMISE-optimal bandwidth for the bias-corrected gamma
#'    kernel, with the unknown density replaced by a gamma reference fitted by
#'    maximum likelihood to the positive observations. The fitted shape is
#'    floored at 3, keeping the fitted mean: the rule needs a shape above 3/2,
#'    and the floor was the most robust choice in simulations with exponential,
#'    lognormal, Weibull and mixture targets. Used for both gamma kernels.
#'
#'    `"RHE"`: A reference rule for symmetric kernels that estimates the
#'    roughness of the density with Hermite polynomials around a fitted
#'    normal density. It uses only the data, not the parametric start.
#'
#'    `"JH"`: Jones & Henderson's normal-reference rule for the Gaussian
#'    copula kernel. Accurate for unimodal data, but it oversmooths
#'    multimodal data badly; capped at 0.5.
#'
#' @section Structure:
#'    The bandwidth selector is a function of four arguments: The data
#'    `x`, a kernel string `kernel`, a start string `start`,
#'    and a support vector `support`. To obtain the functions associated
#'    with these strings, use `get_kernel` and `get_start`. The
#'    function should return a double.
#'
#' @seealso [kdensity()], [stats::bandwidth.kernel()] for the
#'    bandwidth selectors of [stats::density()]. In addition,
#'    [kernels()]; [parametric_starts()]
#'
#' @examples
#' ## Not a serious bandwidth function.
#' silly_width <- function(x, kernel = NULL, start = NULL, support = NULL) {
#'   rexp(1)
#' }
#' kdensity(mtcars$mpg, start = "gumbel", bw = silly_width)
#' @references
#' Jones, M. C., and D. A. Henderson. "Kernel-type density estimation on the unit interval." Biometrika 94.4 (2007): 977-984.
#' Hjort, Nils Lid, and Ingrid K. Glad. "Nonparametric density estimation with a parametric start." The Annals of Statistics (1995): 882-904.
#' Hallberg Szabadváry, Johan. "A Fast, Closed-Form Bandwidth Selector for the Beta Kernel Density Estimator." arXiv preprint arXiv:2601.19553 (2026).
#' Sheather, S. J., and M. C. Jones. "A reliable data-based bandwidth selection method for kernel density estimation." Journal of the Royal Statistical Society: Series B 53.3 (1991): 683-690.
#' Chen, Song Xi. "Probability density function estimation using gamma kernels." Annals of the Institute of Statistical Mathematics 52.3 (2000): 471-480.
#' @name bandwidths
NULL

bw_environment$JH <- function(x, kernel = NULL, start = NULL, support = NULL) {
  ## The data is transfomed through qnorm, with singularities removed.
  transformed_x <- stats::qnorm(x)
  transformed_x <- transformed_x[transformed_x != Inf & transformed_x != -Inf]
  sigma <- stats::sd(transformed_x)
  mu <- mean(transformed_x)
  n <- length(transformed_x)
  min(sigma * (2 * mu^2 * sigma^2 + 3 * (1 - sigma^2)^2)^(-1 / 5) * n^(-1 / 5), 0.5)
}

bw_environment$RHE <- function(x, kernel = NULL, start = NULL, support = NULL) {
  n <- length(x)
  mu <- mean(x)
  sigma <- stats::sd(x)
  z <- (x - mu) / sigma

  ## Probabilist's Hermite polynomials He_j(sqrt(2) * z) for j = 2..5.
  u <- sqrt(2) * z
  u2 <- u * u
  he2 <- u2 - 1
  he3 <- u * (u2 - 3)
  he4 <- u2 * (u2 - 6) + 3
  he5 <- u * (u2 * (u2 - 10) + 15)

  weight <- sqrt(2) * exp(-z^2 / 2)
  delta <- c(
    mean(weight * he2),
    mean(weight * he3),
    mean(weight * he4),
    mean(weight * he5)
  )

  (1 / 4)^(1 / 5) *
    (delta[1]^2 + delta[2]^2 + delta[3]^2 / 2 + delta[4]^2 / 6)^(-1 / 5) *
    sigma * n^(-1 / 5)
}

bw_environment$nrd0 <- function(data, kernel, start, support) stats::bw.nrd0(data)
bw_environment$nrd <- function(data, kernel, start, support) stats::bw.nrd(data)
bw_environment$bcv <- function(data, kernel, start, support) stats::bw.bcv(data)
bw_environment$SJ <- function(data, kernel = NULL, start = NULL, support = NULL) {
  if (!identical(kernel, "gcopula")) {
    return(stats::bw.SJ(data))
  }

  tryCatch(
    compute_probit_sj_bandwidth(data),
    error = function(error) {
      warning(
        paste0(
          "Bandwidth selector 'SJ' failed: ",
          conditionMessage(error),
          ". Falling back to 'JH'."
        ),
        call. = FALSE
      )
      bw_environment$JH(data)
    }
  )
}

bw_environment$HS <- function(x, kernel = NULL, start = NULL, support = NULL) {
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector.")
  }

  if (length(x) < 2L) {
    stop("'x' must have at least 2 observations.")
  }

  if (any(x < 0 | x > 1, na.rm = TRUE)) {
    stop("All values in 'x' must be in [0, 1].")
  }

  tryCatch(
    compute_hs_bandwidth(x),
    error = function(error) {
      warning(
        paste0(
          "Bandwidth selector 'HS' failed: ",
          conditionMessage(error),
          ". Falling back to 'ucv'."
        ),
        call. = FALSE
      )
      bw_environment$ucv(x, kernel = kernel, start = start, support = support)
    }
  )
}

bw_environment$GR <- function(x, kernel = NULL, start = NULL, support = NULL) {
  if (!is.null(kernel) && !kernel %in% c("gamma", "gamma_biased")) {
    stop("Bandwidth selector 'GR' is only defined for the gamma kernels.")
  }

  tryCatch(
    compute_gr_bandwidth(x),
    error = function(error) {
      warning(
        paste0(
          "Bandwidth selector 'GR' failed: ",
          conditionMessage(error),
          ". Falling back to 'ucv'."
        ),
        call. = FALSE
      )
      bw_environment$ucv(x, kernel = kernel, start = start, support = support)
    }
  )
}

bw_environment$ucv <- function(x, kernel = NULL, start = NULL, support = NULL) {
  ## We check for the combination start == "uniform" and kernel == "gaussian",
  ## as this is handled by stats::density's related functions.

  if (!is.null(start)) {
    if (start == "constant" | start == "uniform") {
      if (!is.null(get_kernel(kernel)$sd)) {
        return(stats::bw.ucv(x))
      }
    }
  }


  ## If we are here, we must do our own cross-validation.
  obj_func <- make_ucv_objective(x, kernel, start, support)

  # Use a pilot on the kernel's bandwidth scale. Call the rules directly:
  # their public wrappers can fall back to UCV and would recurse here.
  legacy_pilot <- function() {
    if (kernel %in% c("gcopula", "beta", "beta_biased")) {
      list(value = bw_environment$JH(x), name = "JH")
    } else if (start %in% c("constant", "uniform")) {
      list(value = stats::bw.nrd0(x), name = "nrd0")
    } else {
      list(value = bw_environment$RHE(x), name = "RHE")
    }
  }
  pilot <- tryCatch({
    if (kernel %in% c("gamma", "gamma_biased")) {
      list(value = compute_gr_bandwidth(x), name = "GR")
    } else if (kernel %in% c("beta", "beta_biased")) {
      list(value = suppressWarnings(compute_hs_bandwidth(x)), name = "HS")
    } else if (kernel == "gcopula") {
      list(value = compute_probit_sj_bandwidth(x), name = "SJ")
    } else {
      legacy_pilot()
    }
  }, error = function(e) legacy_pilot())
  using <- pilot$value
  using_str <- pilot$name
  lower <- using / 5
  upper <- 5 * using
  if (kernel %in% c("beta", "beta_biased")) upper <- min(upper, 1 / 4 - 1e-10)
  if (kernel == "gcopula") upper <- min(upper, 1 / 2)

  optimized <- TRUE
  bw <- tryCatch(stats::optimize(obj_func, lower = lower, upper = upper, tol = 0.0001)$minimum,
    error = function(e) {
      optimized <<- FALSE
      warning(paste0("Integration failed when finding bandwidth using 'ucv'. Using '", using_str, "' instead."), call. = FALSE)
      using
    }
  )
  if (optimized && (bw <= lower * 1.01 || bw >= upper * 0.99)) {
    warning(sprintf(
      "'ucv' selected a bandwidth at the edge of its search interval [%g, %g]; the criterion may prefer a bandwidth outside it.",
      lower, upper
    ), call. = FALSE)
  }
  return(bw)
}

## ---------------------------------------------------------------------------
## Internal helpers used by the selectors above.
## ---------------------------------------------------------------------------

## The unbiased cross-validation criterion
##
##   UCV(h) = integral fhat_h^2 - (2 / n) sum_i fhat_{h,-i}(x_i),
##
## estimating ISE(h) - integral f^2. kdensity() normalizes its estimates, so
## both terms use normalized estimates: the leave-one-out sums are divided by
## the full-data normalizer, which approximates each leave-one-out normalizer
## to relative order 1 / n. Leaving the second term unnormalized biases the
## criterion towards whichever bandwidths make the normalizer large; for the
## gamma kernels that is a strong pull towards oversmoothing.
make_ucv_objective <- function(x, kernel, start, support) {
  kernel_obj <- get_kernel(kernel)
  start_obj <- get_start(start)
  kernel_fun <- kernel_obj$kernel
  start_density <- start_obj$density
  start_estimator <- start_obj$estimator
  start_support <- start_obj$support
  n <- length(x)

  # Name of the variable where the density is evaluated. Typically x.
  x_name <- names(formals(start_density))[1]

  full_parameters <- start_estimator(x)

  arguments <- list()
  arguments[[1]] <- x
  names(arguments)[1] <- x_name
  arguments <- append(arguments, as.list(full_parameters))

  dstart <- function(data, parameters) {
    arguments[[1]] <- data
    if (length(parameters) > 0) {
      for (i in seq_along(parameters)) arguments[[i + 1]] <- parameters[[i]]
    }
    do.call(start_density, arguments)
  }

  uniform_start <- start %in% c("constant", "uniform")
  use_compiled_loo <- !uniform_start &&
    !is.null(kernel_obj$compiled_id) &&
    !is.null(start_obj$compiled_id)
  compiled_param_loo <- if (use_compiled_loo) {
    tryCatch(make_compiled_loo_parameters(x, start_obj), error = function(e) NULL)
  } else {
    NULL
  }
  ## Fall back to the R leave-one-out path when the compiled fit fails or
  ## returns non-finite/invalid parameters (degenerate or pathological data),
  ## so a bad fit degrades to the univariateML path and the outer tryCatch
  ## instead of aborting kdensity().
  if (use_compiled_loo &&
      !compiled_loo_parameters_are_valid(
        compiled_param_loo,
        start_obj$compiled_id
      )) {
    use_compiled_loo <- FALSE
  }
  param_loo <- if (uniform_start || use_compiled_loo) {
    NULL
  } else {
    lapply(seq_len(n), function(i) start_estimator(x[-i]))
  }

  parametric_start_vector <- function(data) dstart(data, full_parameters)

  parametric_start_data <- parametric_start_vector(x)
  uniform_data <- rep(1, n)

  function(h) {
    kdensity_sqrd <- kdensity_sq(x,
      h = h,
      kernel_fun = kernel_fun,
      kernel_obj = kernel_obj,
      parametric_start_data = parametric_start_data,
      parametric_start_vector = parametric_start_vector,
      parametric_start = start_density,
      support = support
    )

    term1 <- stats::integrate(kdensity_sqrd,
      lower = support[1],
      upper = support[2]
    )$value

    if (uniform_start && !is.null(kernel_obj$compiled_id)) {
      full_evaluator <- make_kernel_evaluator(
        x = x,
        h = h,
        kernel = kernel_obj,
        parametric_start_data = uniform_data
      )
      full_sum <- sum(full_evaluator(x))
      diagonal <- vapply(x, function(point) {
        kernel_fun(point, point, h) / h
      }, numeric(1))
      term2 <- 2 * (full_sum - mean(diagonal)) / (n - 1)
    } else if (use_compiled_loo) {
      term2 <- compiled_ucv_term2(
        x = x,
        h = h,
        kernel = kernel_obj,
        start = start_obj,
        parameters = compiled_param_loo
      )
    } else {
      term2_vec <- numeric(n)
      for (i in seq_len(n)) {
        # Eval-first: the LOO density at x[i] sets the (possibly asymmetric)
      # kernel's shape from the evaluation point x[i] and evaluates it at the
      # other data points x[-i], consistent with term1 and the estimator.
      term2_vec[i] <- mean(1 / h * kernel_fun(x[i], x[-i], h) /
          dstart(x[-i], param_loo[[i]])) *
          dstart(x[i], param_loo[[i]])
      }
      term2 <- 2 * mean(term2_vec)
    }

    term1 - term2 / attr(kdensity_sqrd, "normalization")
  }
}

kdensity_sq <- function(x, h, kernel_fun, parametric_start, parametric_start_data,
                        parametric_start_vector, support, kernel_obj = NULL) {
  normalization <- 1

  if (is.null(kernel_obj)) {
    kernel_evaluator <- function(y) {
      vapply(y, function(point) {
        mean(1 / h * kernel_fun(point, x, h) / parametric_start_data)
      }, numeric(1))
    }
  } else {
    kernel_evaluator <- make_kernel_evaluator(
      x = x,
      h = h,
      kernel = kernel_obj,
      parametric_start_data = parametric_start_data
    )
  }

  pre_function <- function(y) {
    kernel_evaluator(y) * parametric_start_vector(y)
  }

  normalization <- tryCatch(stats::integrate(pre_function, lower = support[1], upper = support[2])$value,
    error = function(e) {
      stop("Normalization error: The function will not integrate. Two common causes are: 1.) The kernel is non-smooth, try a smooth kernel if possible. 2.) The supplied support is incorrect.")
    }
  )

  return_function <- function(y) {
    parametric_start_vector_y <- parametric_start_vector(y)
    (kernel_evaluator(y) * parametric_start_vector_y / normalization)^2
  }
  attr(return_function, "normalization") <- normalization
  return_function
}

compute_hs_bandwidth <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  x_interior <- x[x > 0 & x < 1]

  if (length(x_interior) == 0L) {
    stop("No data strictly within (0, 1).")
  }

  mu <- mean(x_interior)
  variance <- mean((x_interior - mu)^2)

  if (variance == 0) {
    stop("Sample variance is zero.")
  }

  common <- mu * (1 - mu) / variance - 1
  alpha <- mu * common
  beta <- (1 - mu) * common

  bandwidth <- NA_real_
  use_fallback <- !(alpha > 1.5 && beta > 1.5 && (alpha + beta) > 3)

  if (!use_fallback) {
    log_numerator <- log(2 * alpha + 2 * beta - 5) +
      log(2 * alpha + 2 * beta - 3) +
      lgamma(2 * alpha + 2 * beta - 6) +
      lgamma(alpha) +
      lgamma(beta) +
      lgamma(alpha - 0.5) +
      lgamma(beta - 0.5)

    denominator_term_1 <- (alpha - 1) * (beta - 1)
    denominator_term_2 <- 6 - 4 * beta + alpha * (3 * beta - 4)

    log_denominator <- log(denominator_term_1) +
      log(denominator_term_2) +
      lgamma(2 * alpha - 3) +
      lgamma(2 * beta - 3) +
      lgamma(alpha + beta) +
      lgamma(alpha + beta - 1)

    log_factor <- log(2) + log(n) + 0.5 * log(pi)
    bandwidth <- exp((2 / 5) * (log_numerator - log_denominator - log_factor))
  }

  if (use_fallback) {
    beta_variance <- alpha * beta / ((alpha + beta)^2 * (alpha + beta + 1))
    beta_skewness <- 2 * (beta - alpha) * sqrt(alpha + beta + 1) /
      ((alpha + beta + 2) * sqrt(alpha * beta))
    beta_kurtosis <- 6 * ((alpha - beta)^2 * (alpha + beta + 1) -
      alpha * beta * (alpha + beta + 2)) /
      (alpha * beta * (alpha + beta + 2) * (alpha + beta + 3))

    scale <- sqrt(beta_variance)
    correction <- 1 + abs(beta_skewness) + abs(beta_kurtosis)
    bandwidth <- scale / correction * n^(-0.4)

    warning(
      "MISE rule not applicable; using the HS fallback heuristic.",
      call. = FALSE
    )
  }

  bandwidth
}

compute_gr_bandwidth <- function(x) {
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector.")
  }

  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 2L) {
    stop("'x' must have at least 2 observations.")
  }
  if (any(!is.finite(x)) || any(x < 0)) {
    stop("All values in 'x' must be finite and non-negative.")
  }

  # The gamma reference is fitted to the positive observations; exact zeros
  # have no likelihood under it. The bandwidth still uses the full n.
  x <- x[x > 0]
  if (length(x) < 2L) {
    stop("At least 2 observations must be strictly positive.")
  }

  mean_x <- mean(x)
  mean_log_x <- mean(log(x))
  s <- log(mean_x) - mean_log_x
  if (!is.finite(s) || s <= 0) {
    stop("The gamma reference fit is degenerate.")
  }

  shape <- 1 / (12 * s) *
    (3 - s + sqrt((s - 3)^2 + 24 * s)) * 0.9989 + 0.0102
  tolerance <- .Machine$double.eps^0.25
  for (iteration in seq_len(100L)) {
    step <- (log(shape) - digamma(shape) - s) /
      (1 / shape - trigamma(shape))
    next_shape <- shape - step
    if (!is.finite(next_shape) || next_shape <= 0) {
      stop("The gamma reference fit did not converge.")
    }
    if (abs((shape - next_shape) / shape) < tolerance) {
      shape <- next_shape
      break
    }
    shape <- next_shape
  }
  # The roughness integral below is finite only for shape > 3/2 and blows up
  # as the shape approaches 3/2. Flooring the shape at 3 (keeping the fitted
  # mean) was the most robust choice in the default-bandwidth simulation
  # (tools/default-bandwidth-sim): it tracks GR for gamma-like data and stays
  # sensible for exponential, lognormal and mixture data.
  shape <- max(shape, 3)
  rate <- shape / mean_x

  # For f = Gamma(shape, rate), Chen's two AMISE integrals reduce to
  #
  # V = integral x^(-1/2) f(x) dx
  # R = integral {x f''(x)}^2 dx.
  #
  # Work on the log scale to retain accuracy for concentrated references.
  log_v <- 0.5 * log(rate) + lgamma(shape - 0.5) - lgamma(shape)
  log_r <- 3 * log(rate) + lgamma(2 * shape - 3) -
    (2 * shape - 3) * log(2) - 2 * lgamma(shape) +
    log(shape - 1) + log(3 * shape - 4) - log(4)

  exp((2 / 5) * (
    log_v - log(2 * sqrt(pi)) - log_r - log(n)
  ))
}

## On the probit scale z = qnorm(u), the Gaussian copula estimator is exactly
## a Gaussian kernel estimate of the points rho * qnorm(x_i) with bandwidth
## s = sqrt(1 - rho^2), where rho = 1 - h^2. A Gaussian bandwidth s selected
## for qnorm(x) therefore converts to h = sqrt(1 - sqrt(1 - s^2)), or about
## s / sqrt(2). Observations at 0 or 1 have infinite probits and are dropped,
## as in JH.
compute_probit_sj_bandwidth <- function(x) {
  z <- stats::qnorm(x)
  s <- stats::bw.SJ(z[is.finite(z)])
  s <- min(s, 1)
  sqrt(1 - sqrt(1 - s^2))
}

## ---------------------------------------------------------------------------
## Accessors.
## ---------------------------------------------------------------------------

#' Get bandwidth functions from string.
#'
#' @keywords internal
#' @param bw_str a string specifying the density of interest.
#' @return a bandwidth function.
get_bw <- function(bw_str) {
  assert_(is.character(bw_str))

  bw <- bw_environment[[bw_str]]

  msg <- paste0("The supplied bandwidth function ('", bw_str, "') is not implemented.")
  assert_(!is.null(bw), msg = msg)

  bw
}


#' Get a bandwidth string when 'bw' is unspecified.
#'
#' @keywords internal
#' @param kernel_str a kernel string
#' @param start_str a parametric start string.
#' @param support the support.
#' @return a bandwidth string.

get_standard_bw <- function(kernel_str, start_str, support) {
  ## Every built-in kernel gets a closed-form rule; see the Defaults section
  ## of ?bandwidths and tools/default-bandwidth-sim for the evidence.
  if (kernel_str == "gcopula") {
    return("SJ")
  }
  if (kernel_str %in% c("beta", "beta_biased")) {
    return("HS")
  }
  if (kernel_str %in% c("gamma", "gamma_biased")) {
    return("GR")
  }
  if (is.null(get_kernel(kernel_str)$sd)) {
    return("ucv")
  }
  if (start_str %in% c("constant", "uniform")) "nrd0" else "RHE"
}
