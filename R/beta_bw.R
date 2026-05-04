compute_beta_rot_bandwidth <- function(x) {
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
      "MISE rule not applicable; using the beta_rot fallback heuristic.",
      call. = FALSE
    )
  }

  bandwidth
}
