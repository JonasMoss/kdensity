test_that("FFT grids approximate the supported translation kernels", {
  set.seed(20260720)
  kernels <- c(
    "gaussian", "epanechnikov", "rectangular", "triangular",
    "biweight", "cosine", "optcosine"
  )

  for (kernel_name in kernels) {
    fit <- kdensity(
      rnorm(750),
      kernel = kernel_name,
      start = "uniform",
      bw = 0.25
    )
    grid <- kdensity_grid(
      fit,
      from = -3,
      to = 3,
      n = 257,
      method = "fft",
      fft_n = 4096
    )
    error <- grid$y - fit(grid$x)

    expect_s3_class(grid, "density")
    expect_identical(grid$method, "fft")
    expect_lt(
      sqrt(mean(error^2)),
      if (kernel_name == "rectangular") 2e-3 else 1e-4
    )
  }
})

test_that("FFT grids retain parametric-start corrections", {
  set.seed(20260720)
  fit <- kdensity(
    rnorm(1000, mean = 0.4, sd = 1.3),
    kernel = "gaussian",
    start = "normal",
    bw = 0.2
  )
  grid <- kdensity_grid(
    fit,
    from = -3,
    to = 4,
    n = 301,
    method = "fft",
    fft_n = 4096
  )

  expect_equal(grid$y, fit(grid$x), tolerance = 2e-5)
})

test_that("Gaussian-copula grids use the probit convolution", {
  set.seed(20260720)
  x <- rbeta(1000, 2, 5)

  for (start_name in c("uniform", "beta")) {
    fit <- kdensity(
      x,
      kernel = "gcopula",
      start = start_name,
      bw = 0.12
    )
    grid <- kdensity_grid(
      fit,
      from = 0.005,
      to = 0.995,
      n = 301,
      method = "fft",
      fft_n = 4096
    )

    expect_equal(
      grid$y,
      fit(grid$x),
      tolerance = 5e-5,
      info = start_name
    )
    expect_equal(
      diff(stats::qnorm(range(grid$x))),
      stats::qnorm(0.995) - stats::qnorm(0.005),
      tolerance = 1e-12
    )
  }
})

test_that("grid method selection has an exact fallback", {
  set.seed(20260720)
  fit <- kdensity(rgamma(100, 2), kernel = "gamma", bw = 0.1)

  automatic <- kdensity_grid(fit, n = 101, method = "auto")
  exact <- kdensity_grid(fit, n = 101, method = "exact")

  expect_identical(automatic$method, "exact")
  expect_equal(automatic$x, exact$x)
  expect_equal(automatic$y, exact$y)
  expect_error(
    kdensity_grid(fit, n = 101, method = "fft"),
    "no FFT convolution representation"
  )
})

test_that("Gaussian-copula FFT grids handle support boundaries", {
  set.seed(20260720)
  fit <- kdensity(rbeta(100, 2, 5), kernel = "gcopula", bw = 0.12)

  grid <- kdensity_grid(fit, from = 0, to = 1, method = "fft")
  expect_identical(grid$x[c(1L, length(grid$x))], c(0, 1))

  ## The copula kernel is NaN at the 0/1 endpoints, so pin the finite interior
  ## against exact evaluation rather than comparing NaN == NaN.
  interior <- grid$x > 0 & grid$x < 1
  expect_gt(sum(interior), 2L)
  expect_true(all(is.finite(grid$y[interior])))
  solid <- grid$x > 0.05 & grid$x < 0.95
  expect_equal(grid$y[solid], fit(grid$x[solid]), tolerance = 1e-3)

  expect_error(
    kdensity_grid(fit, from = -0.1, to = 1, method = "fft"),
    "inside \\[0, 1\\]"
  )
})

test_that("grid arguments are validated", {
  fit <- kdensity(precip)

  expect_error(kdensity_grid(precip), "must be a kdensity")
  expect_error(kdensity_grid(fit, n = 1), "at least two")
  expect_error(kdensity_grid(fit, from = 2, to = 1), "from.*to")
  expect_error(kdensity_grid(fit, fft_n = 100), "at least 512")
  expect_error(kdensity_grid(fit, truncate = 0), "positive")
})
