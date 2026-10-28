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

test_that("the gcopula constant does not affect the UCV criterion", {
  # UCV divides the leave-one-out term by the normalizer, so a constant factor
  # in the kernel cancels. Compare the compiled kernel with an R kernel that
  # lacks the constant (the pre-2.0.0 form).
  on.exit(rm(list = "gcopula_unscaled", envir = kernel_environment))
  add_kernel("gcopula_unscaled", list(
    kernel = function(y, x, h) {
      rho <- 1 - h^2
      inside <- rho^2 * (qnorm(y)^2 + qnorm(x)^2) - 2 * rho * qnorm(y) * qnorm(x)
      exp(-inside / (2 * (1 - rho^2)))
    },
    support = c(0, 1)
  ))
  for (start in c("uniform", "beta")) {
    scaled <- make_ucv_objective(gcopula_data, "gcopula", start, c(0, 1))
    unscaled <- make_ucv_objective(gcopula_data, "gcopula_unscaled", start, c(0, 1))
    for (h in c(.12, .25, .4)) {
      expect_equal(scaled(h), unscaled(h), tolerance = 1e-6,
        info = paste(start, h)
      )
    }
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

test_that("gcopula FFT grids agree with exact evaluation", {
  # Includes h > 1 (rho < 0) and data deep in the probit tails, whose terms
  # the FFT path sums exactly.
  extreme_data <- c(gcopula_data, 1e-12, 1e-10, 1 - 1e-12)
  for (data in list(gcopula_data, extreme_data)) {
    for (start in c("uniform", "beta")) {
      for (h in c(.18, .35, 1.2)) {
        fit <- kdensity(data, kernel = "gcopula", start = start, bw = h)
        for (range in list(c(0, 1), c(.05, .95))) {
          grid <- kdensity_grid(fit, from = range[1], to = range[2], n = 65,
                                method = "fft")
          exact <- fit(grid$x)
          info <- paste(length(data), start, h, range[1])
          expect_identical(grid$method, "fft", info = info)
          # A fitted beta start with shape below one is infinite at the
          # boundary, where both paths give 0 * Inf.
          finite <- is.finite(exact)
          expect_identical(is.finite(grid$y), finite, info = info)
          expect_lt(max(abs(grid$y - exact)[finite]) / max(exact[finite]),
                    5e-4)
        }
      }
    }
  }
})

test_that("gcopula FFT grids are exact for the independence copula", {
  # h = 1 gives rho = 0, a constant kernel with no convolution form.
  fit <- kdensity(gcopula_data, kernel = "gcopula", bw = 1)
  grid <- kdensity_grid(fit, from = 0, to = 1, n = 9, method = "fft")
  expect_identical(grid$method, "fft")
  expect_identical(grid$y, fit(grid$x))
})

test_that("gcopula estimates take their limit at 0 and 1", {
  # The probit is infinite at the boundary, where the formula gives Inf - Inf.
  # The kernel vanishes there, except for rho = 0 (h = 1), the constant
  # independence copula.
  for (start in c("uniform", "beta")) {
    fit <- kdensity(gcopula_data, kernel = "gcopula", start = start, bw = .2)
    expect_identical(fit(c(0, 1)), c(0, 0), info = start)
    grid <- kdensity_grid(fit, from = 0, to = 1, n = 33)
    expect_false(anyNA(grid$y), info = start)
    expect_identical(grid$y[c(1L, 33L)], c(0, 0), info = start)
  }

  independence <- kdensity(gcopula_data, kernel = "gcopula", bw = 1)
  expect_equal(independence(c(0, .3, 1)), c(1, 1, 1), tolerance = 1e-12)

  kernel <- get_kernel("gcopula")$kernel
  expect_identical(kernel(0, gcopula_data[1:3], .2), c(0, 0, 0))
  expect_identical(kernel(c(0, 1), .37, .2), c(0, 0))
  expect_identical(kernel(c(0, 1), .37, 1), c(1, 1))
  expect_equal(kernel(.3, .37, .2) / .2, copula_density(.3, .37, .2),
               tolerance = 1e-12)
})
