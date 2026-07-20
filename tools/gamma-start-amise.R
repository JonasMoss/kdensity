# Experimental AMISE selector for the bias-corrected gamma kernel with a
# gamma parametric start.
#
# This file is deliberately outside R/: the selector is research code until
# its boundary expansion and finite-sample behavior have been validated.

laguerre_coefficients <- function(degree, alpha) {
  stopifnot(
    length(degree) == 1L,
    degree >= 0,
    degree == as.integer(degree),
    length(alpha) == 1L,
    alpha > -1
  )

  powers <- 0:degree
  (-1)^powers * exp(
    lgamma(degree + alpha + 1) -
      lgamma(alpha + powers + 1) -
      lgamma(degree - powers + 1) -
      lgamma(powers + 1)
  )
}

polynomial_evaluate <- function(coefficients, x) {
  value <- rep(0, length(x))
  for (coefficient in rev(coefficients)) {
    value <- value * x + coefficient
  }
  value
}

polynomial_add <- function(left, right, right_scale = 1) {
  length_out <- max(length(left), length(right))
  out <- numeric(length_out)
  out[seq_along(left)] <- left
  out[seq_along(right)] <- out[seq_along(right)] + right_scale * right
  out
}

polynomial_product <- function(left, right) {
  out <- numeric(length(left) + length(right) - 1L)
  for (i in seq_along(left)) {
    indices <- i + seq_along(right) - 1L
    out[indices] <- out[indices] + left[[i]] * right
  }
  out
}

gamma_weighted_polynomial_integral <- function(coefficients, power, rate) {
  degrees <- seq_along(coefficients) - 1L
  terms <- coefficients * exp(
    lgamma(power + degrees + 1) -
      (power + degrees + 1) * log(rate)
  )
  sum(terms)
}

gamma_start_fit <- function(x) {
  estimates <- kdensity:::get_start("gamma")$estimator(x)
  c(
    shape = unname(estimates[["shape"]]),
    rate = unname(estimates[["rate"]])
  )
}

gamma_fisher_information <- function(shape, rate) {
  matrix(
    c(
      trigamma(shape), -1 / rate,
      -1 / rate, shape / rate^2
    ),
    nrow = 2L,
    byrow = TRUE
  )
}

gamma_coefficient_score_covariance <- function(coefficients, shape, rate) {
  degrees <- seq_along(coefficients) - 1L
  moments <- exp(lgamma(shape + degrees) - lgamma(shape))
  c(
    sum(coefficients * moments * (digamma(shape + degrees) - digamma(shape))),
    sum(coefficients * moments * (-degrees / rate))
  )
}

gamma_fitted_coefficient_sd <- function(coefficients, shape, rate) {
  score_covariance <- gamma_coefficient_score_covariance(
    coefficients,
    shape,
    rate
  )
  projected_variance <- drop(
    score_covariance %*%
      solve(gamma_fisher_information(shape, rate), score_covariance)
  )
  sqrt(max(1 - projected_variance, 0))
}

#' Experimental Laguerre-reference bandwidth
#'
#' Approximate the density ratio r = f / g in an orthonormal generalized
#' Laguerre basis under a fitted gamma start g. The second derivative of that
#' expansion supplies the leading bias roughness for the bias-corrected gamma
#' kernel. This function is not part of the package API.
#'
#' @param x Strictly positive data.
#' @param degree Largest Laguerre degree retained. Must be at least 2.
#' @param normalized Include the leading effect of kdensity's final
#'   normalization.
#' @param threshold Non-negative multiplier for hard thresholding coefficients
#'   of degree two and above at `threshold / sqrt(n)`. Zero disables
#'   thresholding.
#' @param shrinkage Coefficient shrinkage rule: `"none"`, `"hard"`, `"soft"`,
#'   or positive-part `"garrote"`. When omitted, a positive `threshold`
#'   selects hard thresholding and zero selects none.
#' @param threshold_scale Use the nominal orthonormal-basis standard error or
#'   its parametric-null value after projecting out fitted gamma score
#'   directions.
#' @return A list containing the bandwidth and all intermediate quantities.
laguerre_gamma_bandwidth <- function(x, degree = 4L, normalized = TRUE,
                                     threshold = 0, shrinkage = NULL,
                                     threshold_scale = c(
                                       "orthonormal", "fitted"
                                     )) {
  if (!is.numeric(x) || length(x) < 3L || any(!is.finite(x)) || any(x <= 0)) {
    stop("'x' must contain at least three finite, strictly positive values.")
  }
  degree <- as.integer(degree)
  if (length(degree) != 1L || is.na(degree) || degree < 2L) {
    stop("'degree' must be a single integer of at least 2.")
  }
  if (length(threshold) != 1L || !is.finite(threshold) || threshold < 0) {
    stop("'threshold' must be a single non-negative number.")
  }
  if (is.null(shrinkage)) {
    shrinkage <- if (threshold > 0) "hard" else "none"
  }
  shrinkage <- match.arg(shrinkage, c("none", "hard", "soft", "garrote"))
  threshold_scale <- match.arg(threshold_scale)

  fit <- gamma_start_fit(x)
  shape <- fit[["shape"]]
  rate <- fit[["rate"]]
  if (!is.finite(shape) || !is.finite(rate) || shape <= 0.5 || rate <= 0) {
    stop("The fitted gamma start must have shape > 1/2 and positive rate.")
  }

  transformed_x <- rate * x
  ratio_coefficients <- 1
  ratio_second_coefficients <- 0
  theta <- numeric(degree + 1L)
  theta_used <- numeric(degree + 1L)
  coefficient_sd <- rep(1, degree + 1L)
  coefficient_score_covariance <- matrix(0, nrow = degree + 1L, ncol = 2L)
  theta[[1]] <- 1
  theta_used[[1]] <- 1
  coefficient_sd[[1]] <- 0

  for (order in seq_len(degree)) {
    laguerre <- laguerre_coefficients(order, shape - 1)
    squared_norm <- exp(
      lgamma(order + shape) -
        lgamma(order + 1) -
        lgamma(shape)
    )
    orthonormal <- laguerre / sqrt(squared_norm)
    theta[[order + 1L]] <- mean(
      polynomial_evaluate(orthonormal, transformed_x)
    )
    theta_used[[order + 1L]] <- theta[[order + 1L]]
    if (threshold_scale == "fitted") {
      score_covariance <- gamma_coefficient_score_covariance(
        orthonormal,
        shape,
        rate
      )
      coefficient_score_covariance[order + 1L, ] <- score_covariance
      projected_variance <- drop(
        score_covariance %*%
          solve(gamma_fisher_information(shape, rate), score_covariance)
      )
      coefficient_sd[[order + 1L]] <- sqrt(
        max(1 - projected_variance, 0)
      )
    }
    if (order >= 2L && shrinkage != "none") {
      cutoff <- threshold * coefficient_sd[[order + 1L]] / sqrt(length(x))
      coefficient <- theta_used[[order + 1L]]
      theta_used[[order + 1L]] <- switch(shrinkage,
        hard = if (abs(coefficient) < cutoff) 0 else coefficient,
        soft = sign(coefficient) * max(abs(coefficient) - cutoff, 0),
        garrote = if (coefficient == 0) {
          0
        } else {
          coefficient * max(1 - cutoff^2 / coefficient^2, 0)
        }
      )
    }
    ratio_coefficients <- polynomial_add(
      ratio_coefficients,
      theta_used[[order + 1L]] * orthonormal
    )

    if (order >= 2L) {
      # d^2 L_k^(shape - 1)(t) / dt^2 = L_(k - 2)^(shape + 1)(t).
      second_derivative <- laguerre_coefficients(order - 2L, shape + 1) /
        sqrt(squared_norm)
      ratio_second_coefficients <- polynomial_add(
        ratio_second_coefficients,
        theta_used[[order + 1L]] * second_derivative
      )
    }
  }

  # With t = rate*x,
  #
  #   D(x) = x g(x) r''(x)
  #        = rate^2/Gamma(shape) * exp(-t) t^shape S(t),
  #
  # where S is ratio_second_coefficients. Its integral is rate*a0.
  degrees_second <- seq_along(ratio_second_coefficients) - 1L
  a0 <- sum(
    ratio_second_coefficients * exp(
      lgamma(shape + degrees_second + 1) - lgamma(shape)
    )
  )

  bias_polynomial <- c(0, ratio_second_coefficients)
  if (normalized) {
    bias_polynomial <- polynomial_add(
      bias_polynomial,
      ratio_coefficients,
      right_scale = -a0
    )
  }

  squared_bias_polynomial <- polynomial_product(
    bias_polynomial,
    bias_polynomial
  )
  dimensionless_roughness <- gamma_weighted_polynomial_integral(
    squared_bias_polynomial,
    power = 2 * shape - 2,
    rate = 2
  )
  roughness <- rate^3 / exp(2 * lgamma(shape)) *
    dimensionless_roughness

  # Chen's integrated-variance factor V = integral x^(-1/2) f(x) dx is
  # estimated directly rather than through the truncated ratio expansion.
  variance_factor <- mean(x^(-0.5))
  if (!is.finite(roughness) || roughness <= 0) {
    bandwidth <- Inf
  } else {
    bandwidth <- (
      variance_factor /
        (2 * sqrt(pi) * length(x) * roughness)
    )^(2 / 5)
  }

  coefficient_covariance <- diag(degree + 1L)
  coefficient_covariance[1L, 1L] <- 0
  if (threshold_scale == "fitted") {
    orders <- 2:(degree + 1L)
    score_covariance <- coefficient_score_covariance[orders, , drop = FALSE]
    coefficient_covariance[orders, orders] <- diag(degree) -
      score_covariance %*%
        solve(gamma_fisher_information(shape, rate), t(score_covariance))
    coefficient_covariance <- (
      coefficient_covariance + t(coefficient_covariance)
    ) / 2
  }

  list(
    bandwidth = bandwidth,
    degree = degree,
    normalized = normalized,
    threshold = threshold,
    shrinkage = shrinkage,
    threshold_scale = threshold_scale,
    coefficient_sd = stats::setNames(
      coefficient_sd,
      paste0("theta", 0:degree)
    ),
    coefficient_covariance = coefficient_covariance,
    shape = shape,
    rate = rate,
    theta = stats::setNames(theta, paste0("theta", 0:degree)),
    theta_used = stats::setNames(theta_used, paste0("theta", 0:degree)),
    ratio_coefficients = ratio_coefficients,
    ratio_second_coefficients = ratio_second_coefficients,
    bias_polynomial = bias_polynomial,
    normalization_bias = rate * a0,
    roughness = roughness,
    variance_factor = variance_factor
  )
}
