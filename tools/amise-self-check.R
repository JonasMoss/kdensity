#!/usr/bin/env Rscript

# Deterministic checks for the experimental gamma/beta AMISE helpers.

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1]]))
script_directory <- dirname(script_path)
repository <- dirname(script_directory)

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("The self-check requires the development package 'pkgload'.")
}
pkgload::load_all(repository, quiet = TRUE)
source(file.path(script_directory, "gamma-start-amise.R"))
source(file.path(script_directory, "beta-start-amise.R"))
source(file.path(script_directory, "mise-oracle-helpers.R"))

probabilities <- (seq_len(200L) - 0.5) / 200

gamma_x <- stats::qgamma(probabilities, shape = 3, rate = 1.5)
gamma_result <- laguerre_gamma_bandwidth(
  gamma_x,
  degree = 6L,
  threshold_scale = "fitted"
)
gamma_g <- function(x) {
  stats::dgamma(
    x,
    shape = gamma_result$shape,
    rate = gamma_result$rate
  )
}
gamma_ratio <- function(x) {
  polynomial_evaluate(
    gamma_result$ratio_coefficients,
    gamma_result$rate * x
  )
}
gamma_second <- function(x) {
  gamma_result$rate^2 * polynomial_evaluate(
    gamma_result$ratio_second_coefficients,
    gamma_result$rate * x
  )
}
gamma_a <- stats::integrate(
  function(x) x * gamma_g(x) * gamma_second(x),
  0,
  Inf,
  rel.tol = 1e-10
)$value
gamma_j <- stats::integrate(
  function(x) {
    density <- gamma_g(x) * gamma_ratio(x)
    bias <- x * gamma_g(x) * gamma_second(x) - gamma_a * density
    bias^2
  },
  0,
  Inf,
  rel.tol = 1e-9
)$value

stopifnot(
  isTRUE(all.equal(
    gamma_a,
    gamma_result$normalization_bias,
    tolerance = 1e-8
  )),
  isTRUE(all.equal(
    gamma_j,
    gamma_result$roughness,
    tolerance = 1e-7
  ))
)

beta_x <- stats::qbeta(probabilities, shape1 = 3, shape2 = 5)
beta_result <- jacobi_beta_bandwidth(
  beta_x,
  degree = 6L,
  threshold_scale = "fitted"
)
beta_g <- function(x) {
  stats::dbeta(
    x,
    shape1 = beta_result$shape1,
    shape2 = beta_result$shape2
  )
}
beta_ratio <- function(x) {
  jacobi_polynomial_evaluate(beta_result$ratio_coefficients, x)
}
beta_second <- function(x) {
  jacobi_polynomial_evaluate(
    beta_result$ratio_second_coefficients,
    x
  )
}
beta_a <- stats::integrate(
  function(x) x * (1 - x) * beta_g(x) * beta_second(x),
  0,
  1,
  rel.tol = 1e-10
)$value
beta_j <- stats::integrate(
  function(x) {
    density <- beta_g(x) * beta_ratio(x)
    bias <- x * (1 - x) * beta_g(x) * beta_second(x) - beta_a * density
    bias^2
  },
  0,
  1,
  rel.tol = 1e-9
)$value

stopifnot(
  isTRUE(all.equal(
    beta_a,
    beta_result$normalization_bias,
    tolerance = 1e-8
  )),
  isTRUE(all.equal(
    beta_j,
    beta_result$roughness,
    tolerance = 1e-7
  ))
)

hs_package <- kdensity:::compute_hs_bandwidth(beta_x)
hs_rederived <- normalized_beta_reference_bandwidth(
  beta_x,
  normalized = FALSE
)$bandwidth
stopifnot(isTRUE(all.equal(hs_package, hs_rederived, tolerance = 1e-10)))

gamma_orders <- 2:(gamma_result$degree + 1L)
beta_orders <- 2:(beta_result$degree + 1L)
stopifnot(
  isTRUE(all.equal(
    diag(gamma_result$coefficient_covariance)[gamma_orders],
    unname(gamma_result$coefficient_sd[gamma_orders])^2,
    tolerance = 1e-10
  )),
  isTRUE(all.equal(
    diag(beta_result$coefficient_covariance)[beta_orders],
    unname(beta_result$coefficient_sd[beta_orders])^2,
    tolerance = 1e-10
  ))
)

gamma_test <- apply_global_series_test(gamma_result, length(gamma_x))
beta_test <- apply_global_series_test(beta_result, length(beta_x))
stopifnot(
  is.finite(gamma_test$global_p_value),
  gamma_test$global_p_value >= 0,
  gamma_test$global_p_value <= 1,
  is.finite(beta_test$global_p_value),
  beta_test$global_p_value >= 0,
  beta_test$global_p_value <= 1
)

set.seed(20260720)
saved_seed <- .Random.seed
gamma_bootstrap_test <- apply_bootstrap_global_series_test(
  gamma_result,
  length(gamma_x),
  simulate = function(size) {
    stats::rgamma(
      size,
      shape = gamma_result$shape,
      rate = gamma_result$rate
    )
  },
  refit = function(x) {
    laguerre_gamma_bandwidth(
      x,
      degree = gamma_result$degree,
      threshold_scale = "fitted"
    )
  },
  replicates = 19L
)
stopifnot(
  identical(.Random.seed, saved_seed),
  gamma_bootstrap_test$bootstrap_successful == 19L,
  gamma_bootstrap_test$bootstrap_failures == 0L,
  gamma_bootstrap_test$global_p_value >= 1 / 20,
  gamma_bootstrap_test$global_p_value <= 1
)

beta_bootstrap_test <- apply_bootstrap_global_series_test(
  beta_result,
  length(beta_x),
  simulate = function(size) {
    stats::rbeta(
      size,
      shape1 = beta_result$shape1,
      shape2 = beta_result$shape2
    )
  },
  refit = function(x) {
    jacobi_beta_bandwidth(
      x,
      degree = beta_result$degree,
      threshold_scale = "fitted"
    )
  },
  replicates = 19L
)
stopifnot(
  identical(.Random.seed, saved_seed),
  beta_bootstrap_test$bootstrap_successful == 19L,
  beta_bootstrap_test$bootstrap_failures == 0L,
  beta_bootstrap_test$global_p_value >= 1 / 20,
  beta_bootstrap_test$global_p_value <= 1
)

gamma_smooth_test <- smooth_bootstrap_bandwidth(
  gamma_bootstrap_test,
  power = 1
)
stopifnot(
  gamma_smooth_test$smooth_weight >= 0,
  gamma_smooth_test$smooth_weight <= 1,
  gamma_smooth_test$bandwidth >=
    gamma_smooth_test$unregularized_bandwidth
)

parametric_beta_fit <- kdensity::kdensity(
  beta_x, kernel = "beta", start = "beta", bw = Inf
)
finite_beta_fit <- kdensity::kdensity(
  beta_x, kernel = "beta", start = "beta", bw = 0.1
)
beta_blend_test <- weighted_density_fit(
  parametric_beta_fit,
  finite_beta_fit,
  weight = 0.25
)
evaluation_points <- c(0.2, 0.5, 0.8)
stopifnot(
  isTRUE(all.equal(
    beta_blend_test(evaluation_points),
    0.75 * parametric_beta_fit(evaluation_points) +
      0.25 * finite_beta_fit(evaluation_points),
    tolerance = 1e-12
  )),
  identical(beta_blend_test$bw, finite_beta_fit$bw)
)

message("Experimental AMISE self-checks passed.")
