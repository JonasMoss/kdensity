# Experimental AMISE selector for the bias-corrected beta kernel with a beta
# parametric start.
#
# This file is deliberately outside R/: the selector is research code until
# its boundary expansion and finite-sample behavior have been validated.

jacobi_polynomial_add <- function(left, right, right_scale = 1) {
  length_out <- max(length(left), length(right))
  out <- numeric(length_out)
  out[seq_along(left)] <- left
  out[seq_along(right)] <- out[seq_along(right)] + right_scale * right
  out
}

jacobi_polynomial_product <- function(left, right) {
  out <- numeric(length(left) + length(right) - 1L)
  for (i in seq_along(left)) {
    indices <- i + seq_along(right) - 1L
    out[indices] <- out[indices] + left[[i]] * right
  }
  out
}

jacobi_polynomial_evaluate <- function(coefficients, x) {
  value <- rep(0, length(x))
  for (coefficient in rev(coefficients)) {
    value <- value * x + coefficient
  }
  value
}

generalized_binomial <- function(top, bottom) {
  exp(
    lgamma(top + 1) -
      lgamma(bottom + 1) -
      lgamma(top - bottom + 1)
  )
}

# Coefficients, in increasing powers of x, of
#
#   P_degree^(alpha, beta)(2*x - 1).
#
# The parameter names here follow the standard Jacobi convention, in which
# alpha is the exponent of (1-z) and beta is the exponent of (1+z).
shifted_jacobi_coefficients <- function(degree, alpha, beta) {
  stopifnot(
    length(degree) == 1L,
    degree >= 0,
    degree == as.integer(degree),
    length(alpha) == 1L,
    alpha > -1,
    length(beta) == 1L,
    beta > -1
  )

  out <- numeric(degree + 1L)
  for (power_x in 0:degree) {
    power_x_minus_one <- degree - power_x
    multiplier <- generalized_binomial(degree + alpha, power_x) *
      generalized_binomial(degree + beta, power_x_minus_one)

    expansion <- numeric(degree + 1L)
    added_power <- 0:power_x_minus_one
    expansion[power_x + added_power + 1L] <- multiplier *
      choose(power_x_minus_one, added_power) *
      (-1)^(power_x_minus_one - added_power)
    out <- out + expansion
  }
  out
}

beta_weighted_polynomial_integral <- function(coefficients, shape1, shape2) {
  degrees <- seq_along(coefficients) - 1L
  sum(coefficients * exp(lbeta(shape1 + degrees, shape2)))
}

beta_start_fit <- function(x) {
  estimates <- kdensity:::get_start("beta")$estimator(x)
  c(
    shape1 = unname(estimates[["shape1"]]),
    shape2 = unname(estimates[["shape2"]])
  )
}

beta_fisher_information <- function(shape1, shape2) {
  common <- trigamma(shape1 + shape2)
  matrix(
    c(
      trigamma(shape1) - common, -common,
      -common, trigamma(shape2) - common
    ),
    nrow = 2L,
    byrow = TRUE
  )
}

beta_coefficient_score_covariance <- function(coefficients, shape1, shape2) {
  degrees <- seq_along(coefficients) - 1L
  moments <- exp(
    lbeta(shape1 + degrees, shape2) - lbeta(shape1, shape2)
  )
  score_shape1 <- moments * (
    digamma(shape1 + shape2) -
      digamma(shape1) +
      digamma(shape1 + degrees) -
      digamma(shape1 + shape2 + degrees)
  )
  score_shape2 <- moments * (
    digamma(shape1 + shape2) -
      digamma(shape1 + shape2 + degrees)
  )
  c(
    sum(coefficients * score_shape1),
    sum(coefficients * score_shape2)
  )
}

beta_fitted_coefficient_sd <- function(coefficients, shape1, shape2) {
  score_covariance <- beta_coefficient_score_covariance(
    coefficients,
    shape1,
    shape2
  )
  projected_variance <- drop(
    score_covariance %*%
      solve(beta_fisher_information(shape1, shape2), score_covariance)
  )
  sqrt(max(1 - projected_variance, 0))
}

#' Experimental normalized beta-reference bandwidth
#'
#' Recompute the HS beta reference rule with the leading effect of kdensity's
#' final normalization included in the bias roughness. The same
#' method-of-moments beta fit as the package's HS selector is used so that only
#' the roughness differs.
#'
#' @param x Data in the unit interval.
#' @param normalized Include the leading effect of final normalization.
#' @return A list containing the bandwidth and reference quantities.
normalized_beta_reference_bandwidth <- function(x, normalized = TRUE) {
  if (!is.numeric(x) || length(x) < 2L || any(!is.finite(x)) ||
      any(x < 0 | x > 1)) {
    stop("'x' must contain at least two finite values in [0, 1].")
  }
  x_interior <- x[x > 0 & x < 1]
  if (!length(x_interior)) {
    stop("No data strictly within (0, 1).")
  }

  mu <- mean(x_interior)
  variance <- mean((x_interior - mu)^2)
  if (!is.finite(variance) || variance <= 0) {
    stop("Sample variance is zero.")
  }
  common <- mu * (1 - mu) / variance - 1
  shape1 <- mu * common
  shape2 <- (1 - mu) * common
  if (!is.finite(shape1) || !is.finite(shape2) ||
      shape1 <= 1.5 || shape2 <= 1.5) {
    stop("The fitted beta reference must have both shapes greater than 3/2.")
  }

  # For f = Beta(shape1, shape2),
  #
  # f''(x) = x^(shape1-3) (1-x)^(shape2-3) P2(x) / B(shape1,shape2).
  one_minus_x_squared <- c(1, -2, 1)
  x_one_minus_x <- c(0, 1, -1)
  x_squared <- c(0, 0, 1)
  second_derivative_polynomial <-
    (shape1 - 2) * (shape1 - 1) * one_minus_x_squared
  second_derivative_polynomial <- jacobi_polynomial_add(
    second_derivative_polynomial,
    x_one_minus_x,
    right_scale = -2 * (shape1 - 1) * (shape2 - 1)
  )
  second_derivative_polynomial <- jacobi_polynomial_add(
    second_derivative_polynomial,
    x_squared,
    right_scale = (shape2 - 2) * (shape2 - 1)
  )

  # Since both beta-reference shapes exceed one, f(0) = f(1) = 0 and
  #
  # A = integral x(1-x) f''(x) dx = f(0) + f(1) - 2 = -2.
  normalization_bias <- -2
  bias_polynomial <- second_derivative_polynomial
  if (normalized) {
    bias_polynomial <- jacobi_polynomial_add(
      bias_polynomial,
      x_one_minus_x,
      right_scale = -normalization_bias
    )
  }

  squared_bias_polynomial <- jacobi_polynomial_product(
    bias_polynomial,
    bias_polynomial
  )
  log_beta_start <- lbeta(shape1, shape2)
  roughness <- beta_weighted_polynomial_integral(
    squared_bias_polynomial,
    shape1 = 2 * shape1 - 3,
    shape2 = 2 * shape2 - 3
  ) / exp(2 * log_beta_start)
  variance_factor <- exp(
    lbeta(shape1 - 0.5, shape2 - 0.5) - log_beta_start
  )
  bandwidth <- (
    variance_factor /
      (2 * sqrt(pi) * length(x) * roughness)
  )^(2 / 5)

  list(
    bandwidth = bandwidth,
    normalized = normalized,
    shape1 = shape1,
    shape2 = shape2,
    bias_polynomial = bias_polynomial,
    normalization_bias = normalization_bias,
    roughness = roughness,
    variance_factor = variance_factor
  )
}

#' Experimental Jacobi-reference bandwidth
#'
#' Approximate the density ratio r = f / g in an orthonormal shifted-Jacobi
#' basis under a fitted beta start g. The second derivative of that expansion
#' supplies the leading bias roughness for the bias-corrected beta kernel.
#' This function is not part of the package API.
#'
#' @param x Data strictly inside the unit interval.
#' @param degree Largest Jacobi degree retained. Must be at least 2.
#' @param normalized Include the leading effect of kdensity's final
#'   normalization.
#' @param threshold Non-negative multiplier for hard thresholding coefficients
#'   of degree two and above at `threshold / sqrt(n)`. Zero disables
#'   thresholding.
#' @param shrinkage Coefficient shrinkage rule: `"none"`, `"hard"`, `"soft"`,
#'   or positive-part `"garrote"`. When omitted, a positive `threshold`
#'   selects hard thresholding and zero selects none.
#' @param threshold_scale Use the nominal orthonormal-basis standard error or
#'   its parametric-null value after projecting out fitted beta score
#'   directions.
#' @return A list containing the bandwidth and all intermediate quantities.
jacobi_beta_bandwidth <- function(x, degree = 4L, normalized = TRUE,
                                  threshold = 0, shrinkage = NULL,
                                  threshold_scale = c(
                                    "orthonormal", "fitted"
                                  )) {
  if (!is.numeric(x) || length(x) < 3L || any(!is.finite(x)) ||
      any(x <= 0 | x >= 1)) {
    stop(
      "'x' must contain at least three finite values strictly inside (0, 1)."
    )
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

  fit <- beta_start_fit(x)
  shape1 <- fit[["shape1"]]
  shape2 <- fit[["shape2"]]
  if (!is.finite(shape1) || !is.finite(shape2) ||
      shape1 <= 0.5 || shape2 <= 0.5) {
    stop("The fitted beta start must have both shapes greater than 1/2.")
  }

  ratio_coefficients <- 1
  ratio_second_coefficients <- 0
  theta <- numeric(degree + 1L)
  theta_used <- numeric(degree + 1L)
  coefficient_sd <- rep(1, degree + 1L)
  coefficient_score_covariance <- matrix(0, nrow = degree + 1L, ncol = 2L)
  theta[[1]] <- 1
  theta_used[[1]] <- 1
  coefficient_sd[[1]] <- 0

  # A Beta(shape1, shape2) density corresponds to the shifted Jacobi weight
  # with standard parameters alpha = shape2 - 1 and beta = shape1 - 1.
  for (order in seq_len(degree)) {
    jacobi <- shifted_jacobi_coefficients(
      order,
      alpha = shape2 - 1,
      beta = shape1 - 1
    )
    log_squared_norm <- lgamma(order + shape1) +
      lgamma(order + shape2) +
      lgamma(shape1 + shape2) -
      log(2 * order + shape1 + shape2 - 1) -
      lgamma(order + 1) -
      lgamma(order + shape1 + shape2 - 1) -
      lgamma(shape1) -
      lgamma(shape2)
    orthonormal <- jacobi / exp(log_squared_norm / 2)

    theta[[order + 1L]] <- mean(
      jacobi_polynomial_evaluate(orthonormal, x)
    )
    theta_used[[order + 1L]] <- theta[[order + 1L]]
    if (threshold_scale == "fitted") {
      score_covariance <- beta_coefficient_score_covariance(
        orthonormal,
        shape1,
        shape2
      )
      coefficient_score_covariance[order + 1L, ] <- score_covariance
      projected_variance <- drop(
        score_covariance %*%
          solve(
            beta_fisher_information(shape1, shape2),
            score_covariance
          )
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
    ratio_coefficients <- jacobi_polynomial_add(
      ratio_coefficients,
      theta_used[[order + 1L]] * orthonormal
    )

    if (order >= 2L) {
      # d^2 P_k^(alpha,beta)(2*x-1) / dx^2
      #   = (k+alpha+beta+1)(k+alpha+beta+2)
      #       P_(k-2)^(alpha+2,beta+2)(2*x-1).
      second_derivative <- (
        (order + shape1 + shape2 - 1) *
          (order + shape1 + shape2)
      ) * shifted_jacobi_coefficients(
        order - 2L,
        alpha = shape2 + 1,
        beta = shape1 + 1
      ) / exp(log_squared_norm / 2)
      ratio_second_coefficients <- jacobi_polynomial_add(
        ratio_second_coefficients,
        theta_used[[order + 1L]] * second_derivative
      )
    }
  }

  # D(x) = x(1-x) g(x) r''(x). Its integral A is the leading raw mass
  # error multiplied by 2/b.
  log_beta_start <- lbeta(shape1, shape2)
  normalization_bias <- beta_weighted_polynomial_integral(
    ratio_second_coefficients,
    shape1 = shape1 + 1,
    shape2 = shape2 + 1
  ) / exp(log_beta_start)

  # Factor D - A*f as
  #
  #   x^(shape1-1) (1-x)^(shape2-1) / B(shape1,shape2)
  #     * [x(1-x) S(x) - A R(x)].
  weighted_second <- jacobi_polynomial_product(
    c(0, 1, -1),
    ratio_second_coefficients
  )
  bias_polynomial <- weighted_second
  if (normalized) {
    bias_polynomial <- jacobi_polynomial_add(
      bias_polynomial,
      ratio_coefficients,
      right_scale = -normalization_bias
    )
  }

  squared_bias_polynomial <- jacobi_polynomial_product(
    bias_polynomial,
    bias_polynomial
  )
  roughness <- beta_weighted_polynomial_integral(
    squared_bias_polynomial,
    shape1 = 2 * shape1 - 1,
    shape2 = 2 * shape2 - 1
  ) / exp(2 * log_beta_start)

  # Chen's integrated-variance factor
  # V = integral [x(1-x)]^(-1/2) f(x) dx is estimated directly rather than
  # through the truncated density-ratio expansion.
  variance_factor <- mean((x * (1 - x))^(-0.5))
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
        solve(
          beta_fisher_information(shape1, shape2),
          t(score_covariance)
        )
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
    shape1 = shape1,
    shape2 = shape2,
    theta = stats::setNames(theta, paste0("theta", 0:degree)),
    theta_used = stats::setNames(theta_used, paste0("theta", 0:degree)),
    ratio_coefficients = ratio_coefficients,
    ratio_second_coefficients = ratio_second_coefficients,
    bias_polynomial = bias_polynomial,
    normalization_bias = normalization_bias,
    roughness = roughness,
    variance_factor = variance_factor
  )
}
