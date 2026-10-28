gcopula_data <- c(.08, .13, .21, .29, .37, .48, .56, .65, .73, .84, .92)

copula_density <- function(y, x, h) {
  rho <- 1 - h^2
  sigma <- sqrt(1 - rho^2)
  stats::dnorm((stats::qnorm(y) - rho * stats::qnorm(x)) / sigma) /
    (sigma * stats::dnorm(stats::qnorm(y)))
}

test_that("gcopula kernels include the copula density constant", {
  y <- c(.1, .3, .5, .7, .9)
  for (h in c(.18, .35)) {
    expected <- vapply(y, function(z) {
      mean(copula_density(z, gcopula_data, h))
    }, numeric(1))
    fit <- kdensity(gcopula_data, kernel = "gcopula", bw = h,
                    normalized = FALSE)
    expect_equal(fit(y), expected, tolerance = 1e-12)
    expect_equal(stats::integrate(fit, 0, 1)$value, 1, tolerance = 1e-8)
    expect_equal(get_kernel("gcopula")$kernel(y, .37, h) / h,
                 copula_density(y, .37, h), tolerance = 1e-12)
    expect_equal(kdensity(gcopula_data, kernel = "gcopula", bw = h)(y),
                 expected, tolerance = 1e-12)
  }
})

test_that("gcopula normalized beta starts agree with the density formula", {
  for (h in c(.18, .35)) {
    fit <- kdensity(gcopula_data, kernel = "gcopula", start = "beta", bw = h)
    start <- environment(fit)$parametric_start_vector
    raw <- function(y) vapply(y, function(z) {
      start(z) * mean(copula_density(z, gcopula_data, h) / start(gcopula_data))
    }, numeric(1))
    normalizer <- stats::integrate(raw, 0, 1)$value
    y <- c(.1, .3, .5, .7, .9)
    expect_equal(fit(y), raw(y) / normalizer, tolerance = 1e-12)
    grid <- kdensity_grid(fit, from = .05, to = .95, n = 33, method = "exact")
    expect_equal(grid$y, raw(grid$x) / normalizer, tolerance = 1e-12)
  }
})

test_that("gcopula UCV retains the reference bandwidth", {
  for (start in c("uniform", "beta")) {
    fit <- kdensity(gcopula_data, kernel = "gcopula", start = start, bw = "ucv")
    expect_equal(fit$h, 0.2499433619183222, tolerance = 1e-12)
  }
})

test_that("gcopula rejects boundary observations before fitting", {
  for (boundary in c(0, 1)) {
    expect_error(kdensity(c(boundary, gcopula_data), kernel = "gcopula"),
                 "strictly between 0 and 1.*beta kernel.*off the boundary")
  }
  expect_error(kdensity(c(0, gcopula_data), support = c(0, 1)),
               "strictly between 0 and 1")
})

test_that("gcopula FFT grids retain reference values", {
  # Central points from the unmodified package's 33-point probit grids.
  reference <- c(1.00358810943314, 1.02794903884636,
                 1.01901128951024, 1.08449363545763)
  i <- 0L
  for (start in c("uniform", "beta")) {
    for (h in c(.18, .35)) {
      i <- i + 1L
      fit <- kdensity(gcopula_data, kernel = "gcopula", start = start, bw = h)
      grid <- kdensity_grid(fit, n = 33, method = "fft")
      expect_equal(grid$y[17], reference[i], tolerance = 1e-12)
    }
  }
})
