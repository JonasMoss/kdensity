expect_error(get_bw(5))
expect_error(get_bw("no_kernel"))
expect_silent(get_bw("HS"))
expect_silent(get_bw("GR"))

test_that("every built-in kernel defaults to a closed-form rule", {
  expect_equal(get_standard_bw("gcopula", "uniform"), "SJ")
  expect_equal(get_standard_bw("gcopula", "constant"), "SJ")
  expect_equal(get_standard_bw("gcopula", "beta"), "SJ")
  expect_equal(get_standard_bw("beta", "uniform"), "HS")
  expect_equal(get_standard_bw("beta", "constant"), "HS")
  expect_equal(get_standard_bw("beta", "gaussian"), "HS")
  expect_equal(get_standard_bw("beta_biased", "beta"), "HS")
  expect_equal(get_standard_bw("gamma", "normal"), "GR")
  expect_equal(get_standard_bw("gamma", "uniform"), "GR")
  expect_equal(get_standard_bw("gamma_biased", "gamma"), "GR")
  expect_equal(get_standard_bw("normal", "normal"), "RHE")
  expect_equal(get_standard_bw("uniform", "normal"), "RHE")
  expect_equal(get_standard_bw("epanechnikov", "constant"), "nrd0")
  expect_equal(get_standard_bw("triangular", "uniform"), "nrd0")
  add_kernel("asymmetric_custom", list(
    kernel = function(y, x, h) stats::dgamma(x, y / h + 1, scale = h) * h,
    support = c(0, Inf)
  ))
  expect_equal(get_standard_bw("asymmetric_custom", "uniform"), "ucv")
})

set.seed(123)
gamma_ucv_sample <- rgamma(60, 2, 2)
## Regression pin for Chen's modified gamma kernel with the normalized UCV
## criterion (2.0.0). For this sample the ISE-optimal bandwidth is about 0.11.
expect_equal(
  get_bw("ucv")(gamma_ucv_sample, "gamma", "uniform", c(0, Inf)),
  0.1077853552,
  tolerance = 1e-8
)

set.seed(1)
regular_beta_sample <- rbeta(200, 2, 5)
hs_value <- get_bw("HS")(regular_beta_sample, "beta", "uniform", c(0, 1))
expect_gt(hs_value, 0)
## Regression pin for the closed-form HS path; recomputed manually from
## the alpha/beta moment fit at seed = 1.
expect_equal(hs_value, 0.0209275823, tolerance = 1e-8)

set.seed(1)
rhe_sample <- rnorm(500)
## Regression pin for the inlined Hermite polynomial RHE selector;
## matches the previous EQL-based implementation to ~1e-12.
expect_equal(get_bw("RHE")(rhe_sample), 0.4314342622, tolerance = 1e-8)

set.seed(1)
boundary_beta_sample <- rbeta(200, 0.5, 0.5)
expect_warning(
  get_bw("HS")(boundary_beta_sample, "beta", "uniform", c(0, 1)),
  "fallback heuristic"
)

expect_error(get_bw("HS")("not numeric", "beta", "uniform", c(0, 1)), "'x' must be a numeric vector.")
expect_error(get_bw("HS")(0.5, "beta", "uniform", c(0, 1)), "at least 2 observations")
expect_error(get_bw("HS")(c(-0.1, 0.2), "beta", "uniform", c(0, 1)), "must be in \\[0, 1\\]")
expect_error(compute_hs_bandwidth(c(0, 1, NA)), "No data strictly within")
expect_error(compute_hs_bandwidth(c(0.4, 0.4, 0.4)), "Sample variance is zero.")

set.seed(42)
gamma_reference_sample <- rgamma(300, shape = 3.2, rate = 1.7)
## Regression pin. The fitted shape (2.84) is below 3, so the floor applies;
## before 2.0.0 this sample gave 0.070709052899177.
expect_equal(
  get_bw("GR")(
    gamma_reference_sample, "gamma", "uniform", c(0, Inf)
  ),
  0.067252214229985,
  tolerance = 1e-10
)
expect_error(
  get_bw("GR")(gamma_reference_sample, "beta", "uniform", c(0, 1)),
  "only defined for the gamma kernels"
)
expect_equal(
  get_bw("GR")(gamma_reference_sample, "gamma_biased", "uniform", c(0, Inf)),
  get_bw("GR")(gamma_reference_sample, "gamma", "uniform", c(0, Inf))
)
expect_error(
  compute_gr_bandwidth(c(-1, 1, 2)),
  "non-negative"
)
expect_error(
  compute_gr_bandwidth(c(0, 0, 1)),
  "At least 2 observations must be strictly positive"
)
expect_error(
  compute_gr_bandwidth("not numeric"),
  "numeric vector"
)
expect_error(
  compute_gr_bandwidth(1),
  "at least 2 observations"
)
expect_error(
  compute_gr_bandwidth(rep(2, 3)),
  "gamma reference fit is degenerate"
)
expect_error(
  compute_gr_bandwidth(c(.Machine$double.xmin, 1)),
  "did not converge"
)
test_that("GR floors the fitted gamma shape at 3 and keeps the mean", {
  # Chen's AMISE bandwidth for a Gamma(shape, rate) reference, written out.
  gr_formula <- function(shape, rate, n) {
    v <- sqrt(rate) * gamma(shape - 0.5) / gamma(shape)
    r <- rate^3 * gamma(2 * shape - 3) / 2^(2 * shape - 3) / gamma(shape)^2 *
      (shape - 1) * (3 * shape - 4) / 4
    (v / (2 * sqrt(pi) * r * n))^(2 / 5)
  }
  # Exponential-like data: fitted shape near 1, floored at 3.
  x <- qgamma(ppoints(100), shape = 1, rate = 2)
  expect_lt(get_start("gamma")$estimator(x)[["shape"]], 3)
  expect_equal(compute_gr_bandwidth(x), gr_formula(3, 3 / mean(x), 100),
    tolerance = 1e-10
  )
  # Above the floor the fitted shape is used as is.
  y <- qgamma(ppoints(200), shape = 6, rate = 2)
  fit <- get_start("gamma")$estimator(y)
  expect_equal(compute_gr_bandwidth(y),
    gr_formula(fit[["shape"]], fit[["shape"]] / mean(y), 200),
    tolerance = 1e-6
  )
  # Exact zeros are dropped from the reference fit but counted in n.
  expect_equal(compute_gr_bandwidth(c(0, y)),
    gr_formula(fit[["shape"]], fit[["shape"]] / mean(y), 201),
    tolerance = 1e-6
  )
})

test_that("SJ for the Gaussian copula kernel works on the probit scale", {
  set.seed(3)
  x <- rbeta(300, 2, 5)
  h <- get_bw("SJ")(x, "gcopula", "uniform", c(0, 1))
  s <- stats::bw.SJ(stats::qnorm(x))
  # The copula kernel's probit-scale sd sqrt(1 - rho^2) equals the SJ bandwidth.
  rho <- 1 - h^2
  expect_equal(sqrt(1 - rho^2), s, tolerance = 1e-12)
  # Observations at the boundary are dropped, as in JH.
  expect_equal(get_bw("SJ")(c(0, x, 1), "gcopula", "uniform", c(0, 1)), h)
  # Other kernels keep stats::bw.SJ on the data scale.
  expect_equal(get_bw("SJ")(x, "gaussian", "uniform", c(-Inf, Inf)), stats::bw.SJ(x))
  # Degenerate data fall back to JH with a warning.
  expect_warning(
    fallback <- get_bw("SJ")(rep(0.4, 5), "gcopula", "uniform", c(0, 1)),
    "Falling back to 'JH'"
  )
  expect_equal(fallback, bw_environment$JH(rep(0.4, 5)))
  # Very wide probit-scale bandwidths are capped at the independence copula.
  expect_equal(compute_probit_sj_bandwidth(pnorm(c(-8, 8))), 1)
})

expect_error(
  kdensity_sq(
    x = c(0.2, 0.8),
    h = 0.1,
    kernel_fun = function(y, x, h) stop("boom"),
    parametric_start = function(x) rep(1, length(x)),
    parametric_start_data = c(1, 1),
    parametric_start_vector = function(y) rep(1, length(y)),
    support = c(0, 1)
  ),
  "Normalization error: The function will not integrate"
)

expect_error(
  kdensity(mtcars$mpg, kernel = "uniform", start = "gamma"),
  "Normalization error: The function will not integrate.Two common causes are: 1.) The kernel is non-smooth, try a smooth kernel if possible. 2.) The supplied support is incorrect."
)

test_that("the UCV objective normalizes both terms", {
  # Independent reference: integral of the squared normalized estimate minus
  # twice the mean leave-one-out estimate, divided by the full-data normalizer.
  ucv_reference <- function(x, h, kernel_name, start_name) {
    kernel <- get_kernel(kernel_name)$kernel
    start <- get_start(start_name)
    dstart <- function(data, estimates) {
      do.call(start$density, c(list(data), as.list(estimates)))
    }
    full <- start$estimator(x)
    unnormalized <- function(y) {
      vapply(y, function(point) {
        mean(kernel(point, x, h) / h / dstart(x, full)) * dstart(point, full)
      }, numeric(1))
    }
    normalization <- integrate(unnormalized, 0, Inf)$value
    term1 <- integrate(function(y) (unnormalized(y) / normalization)^2, 0, Inf)$value
    loo <- vapply(seq_along(x), function(i) {
      estimates <- start$estimator(x[-i])
      mean(kernel(x[i], x[-i], h) / h / dstart(x[-i], estimates)) *
        dstart(x[i], estimates)
    }, numeric(1))
    term1 - 2 * mean(loo) / normalization
  }

  set.seed(20261028)
  x <- rgamma(40, 2, 2)
  cases <- list(
    c("gamma", "uniform"),
    c("gamma_biased", "uniform"),
    c("gamma", "gamma"),
    c("gamma", "invgauss")
  )
  for (case in cases) {
    objective <- make_ucv_objective(x, case[[1]], case[[2]], c(0, Inf))
    for (h in c(0.05, 0.2, 0.6)) {
      expect_equal(objective(h), ucv_reference(x, h, case[[1]], case[[2]]),
        tolerance = 1e-5, info = paste(case[[1]], case[[2]], h)
      )
    }
  }

  # A failing compiled leave-one-out fit must fall back to the R path, which
  # gives the same criterion.
  local_mocked_bindings(
    make_compiled_loo_parameters = function(...) stop("compiled fit failed")
  )
  objective <- make_ucv_objective(x, "gamma", "gamma", c(0, Inf))
  expect_equal(objective(0.2), ucv_reference(x, 0.2, "gamma", "gamma"),
    tolerance = 1e-5
  )
})
