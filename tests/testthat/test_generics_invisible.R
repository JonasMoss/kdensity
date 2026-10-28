obj <- kdensity(precip, start = "gaussian")
expect_identical(print(obj), obj)
expect_identical(summary(obj), obj)
expect_identical(plot(obj), obj)
expect_identical(plot(obj, grid = "exact"), obj)
expect_identical(plot(obj, plot_start = TRUE), obj)
expect_identical(points(obj), obj)
expect_identical(lines(obj), obj)

expect_identical(plot(obj, grid = "auto", zero_line = FALSE), obj)
expect_identical(
  plot(
    obj,
    range = range(obj$x),
    grid = "auto",
    zero_line = FALSE
  ),
  obj
)

set.seed(20282860)
copula <- kdensity(
  rbeta(100, 0.7, 0.7),
  kernel = "gcopula",
  start = "beta",
  bw = 0.04
)
expect_identical(
  plot(
    copula,
    range = seq(0, 1, length.out = 1000),
    zero_line = FALSE
  ),
  copula
)
