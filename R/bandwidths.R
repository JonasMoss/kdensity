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
#' @section Bandwidth selectors:
#'    `"nrd0", "nrd", "bcv", "SJ"`: Bandwidth selectors from `stats`.
#'    They are documented in `[bandwidth][stats::bandwidth] stats:bandwidth`.
#'    "nrd0" is the standard bandwidth selector for symmetric kernels with
#'    constant parametric starts.
#'
#'    `"ucv"`: Unbiased cross validation. The standard option for
#'    asymmetric kernels.
#'
#'    `"HS"`: Hallberg Szabadváry's closed-form reference rule for the
#'    beta kernel on the unit interval. This is the default for
#'    `kernel = "beta"` with a uniform or constant start.
#'
#'    `"RHE"`: Selector for parametric starts with a symmetric kernel,
#'    based on a reference rule with Hermite polynomials.
#'    Described in Hjort & Glad (1995). The default method in `kdensity` when a parametric
#'    start is supplied and the kernel is symmetric.
#'
#'    `"JH"`: Selector for the Gaussian copula kernel, based on
#'    normal reference rule. Described in Jones & Henderson. The default method when
#'    the `gcopula` kernel is used in `kdensity`.
#'
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
bw_environment$SJ <- function(data, kernel, start, support) stats::bw.SJ(data)

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
      for (i in 1:length(parameters)) arguments[[i + 1]] <- parameters[[i]]
    }
    do.call(start_density, arguments)
  }

  param_loo <- lapply(1:n, function(i) start_estimator(x[-i]))

  parametric_start_vector <- function(data) dstart(data, full_parameters)

  parametric_start_data <- parametric_start_vector(x)

  obj_func <- function(h) {
    kdensity_sqrd <- kdensity_sq(x,
      h = h,
      kernel_fun = kernel_fun,
      parametric_start_data = parametric_start_data,
      parametric_start_vector = parametric_start_vector,
      parametric_start = start_density,
      support = support
    )

    term1 <- stats::integrate(kdensity_sqrd,
      lower = support[1],
      upper = support[2]
    )$value

    term2_vec <- rep(0, n)
    for (i in 1:n) {
      # Will not work fro asymmetric? Difference between x and y in kernel.
      term2_vec[i] <- mean(1 / h * kernel_fun(x[-i], x[i], h) /
        dstart(x[-i], param_loo[[i]])) *
        dstart(x[i], param_loo[[i]])
    }
    term2 <- 2 * mean(term2_vec)

    obj_func_value <- term1 - term2
    return(obj_func_value)
  }


  # The range of allowable bandwidths vary from kernel to kernel.

  eps <- 10^-10

  if (kernel == "gcopula" | kernel == "beta" | kernel == "beta_biased") {
    using_str <- "JH"
    using <- bw_environment$JH(x)
    lower <- 1 / 4 * using
    upper <- 1 / 4 - eps
  } else if (start == "constant" | start == "uniform") {
    using_str <- "nrd0"
    using <- stats::bw.nrd0(x)
    lower <- 1 / 5 * using
    upper <- 5 * using
  } else {
    using_str <- "RHE"
    using <- bw_environment$RHE(x)
    lower <- 1 / 5 * using
    upper <- 5 * using
  }

  bw <- tryCatch(stats::optimize(obj_func, lower = lower, upper = upper, tol = 0.0001)$minimum,
    error = function(e) {
      warning(paste0("Integration failed when finding bandwidth using 'ucv'. Using '", using_str, "' instead."), call. = FALSE)
      using
    }
  )
  return(bw)
}

## ---------------------------------------------------------------------------
## Internal helpers used by the selectors above.
## ---------------------------------------------------------------------------

kdensity_sq <- function(x, h, kernel_fun, parametric_start, parametric_start_data,
                        parametric_start_vector, support) {
  normalization <- 1

  pre_function <- function(y) {
    sapply(y, function(y) mean(1 / h * kernel_fun(y, x, h) / parametric_start_data) * parametric_start_vector(y))
  }

  normalization <- tryCatch(stats::integrate(pre_function, lower = support[1], upper = support[2])$value,
    error = function(e) {
      stop("Normalization error: The function will not integrate. Two common causes are: 1.) The kernel is non-smooth, try a smooth kernel if possible. 2.) The supplied support is incorrect.")
    }
  )

  return_function <- function(y) {
    n <- length(y)
    parametric_start_vector_y <- parametric_start_vector(y)
    sapply(1:n, function(i) {
      (1 / h * mean(kernel_fun(y[i], x, h) * parametric_start_vector_y[i] / parametric_start_data) / normalization)^2
    })
  }
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
  if (kernel_str == "gcopula" & (start_str == "constant" |
    start_str == "uniform")) {
    bw <- "JH"
  } else if (kernel_str == "beta" & (start_str == "constant" |
    start_str == "uniform")) {
    bw <- "HS"
  } else if (start_str != "constant" & start_str != "uniform") {
    if (!is.null(get_kernel(kernel_str)$sd)) {
      bw <- "RHE"
    } else {
      bw <- "ucv"
    }
  } else if (!is.null(get_kernel(kernel_str)$sd)) {
    bw <- "nrd0"
  } else {
    bw <- "ucv"
  }

  bw
}
