context("bandwidths")

t <- function(x) {}
good <- function(x, kernel, start, support) {}

expect_error(add_bw("t", t))
expect_silent(add_bw("t", good))
expect_error(get_bw(5))
expect_error(get_bw("no_kernel"))

expect_equal(get_standard_bw(kernel_str = "gcopula", start_str = "uniform"), "JH")
expect_equal(get_standard_bw(kernel_str = "gcopula", start_str = "constant"), "JH")
expect_equal(get_standard_bw(kernel_str = "beta", start_str = "constant"), "ucv")
expect_equal(get_standard_bw(kernel_str = "beta", start_str = "gaussian"), "ucv")
expect_equal(get_standard_bw(kernel_str = "normal", start_str = "normal"), "RHE")
expect_equal(get_standard_bw(kernel_str = "uniform", start_str = "normal"), "RHE")
expect_equal(get_standard_bw(kernel_str = "gamma", start_str = "normal"), "ucv")
expect_equal(get_standard_bw(kernel_str = "gamma", start_str = "uniform"), "ucv")
expect_equal(get_standard_bw(kernel_str = "epanechnikov", start_str = "constant"), "nrd0")
expect_equal(get_standard_bw(kernel_str = "triangular", start_str = "uniform"), "nrd0")

expect_error(
  kdensity(mtcars$mpg, kernel = "uniform", start = "gamma"),
  "Normalization error: The function will not integrate.Two common causes are: 1.) The kernel is non-smooth, try a smooth kernel if possible. 2.) The supplied support is incorrect."
)
