context("bandwidths")

expect_error(get_bw(5))
expect_error(get_bw("no_kernel"))
expect_silent(get_bw("HS"))

expect_equal(get_standard_bw(kernel_str = "gcopula", start_str = "uniform"), "JH")
expect_equal(get_standard_bw(kernel_str = "gcopula", start_str = "constant"), "JH")
expect_equal(get_standard_bw(kernel_str = "beta", start_str = "uniform"), "HS")
expect_equal(get_standard_bw(kernel_str = "beta", start_str = "constant"), "HS")
expect_equal(get_standard_bw(kernel_str = "beta", start_str = "gaussian"), "ucv")
expect_equal(get_standard_bw(kernel_str = "normal", start_str = "normal"), "RHE")
expect_equal(get_standard_bw(kernel_str = "uniform", start_str = "normal"), "RHE")
expect_equal(get_standard_bw(kernel_str = "gamma", start_str = "normal"), "ucv")
expect_equal(get_standard_bw(kernel_str = "gamma", start_str = "uniform"), "ucv")
expect_equal(get_standard_bw(kernel_str = "epanechnikov", start_str = "constant"), "nrd0")
expect_equal(get_standard_bw(kernel_str = "triangular", start_str = "uniform"), "nrd0")

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
