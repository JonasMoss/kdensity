# Degenerate / graceful-failure paths. The compiled leave-one-out fits report
# non-convergence the "Result<>" way (an empty std::optional -> NaN) instead of
# throwing, and the R layer falls back rather than aborting.

test_that("compiled leave-one-out fits report degeneracy with NaN, not errors", {
  ## n < 3: too few retained observations for a leave-one-out fit.
  m <- designated_loo_parameters_cpp(c(1, 2), 4L)
  expect_equal(dim(m), c(2L, 2L))
  expect_true(all(is.na(m)))

  ## Constant beta data: the geometric means make the initial guess non-finite.
  mb <- designated_loo_parameters_cpp(rep(0.5, 5), 6L)
  expect_true(all(is.na(mb)))

  ## Extreme gamma data: the Newton step hits NaN for some leave-one-out subsets.
  mg <- designated_loo_parameters_cpp(c(1e-300, 1, 2, 3, 4), 4L)
  expect_true(anyNA(mg))

  ## The structural assertion (a programming error, unreachable from R) is the
  ## one remaining exception.
  expect_error(designated_loo_parameters_cpp(c(1, 2, 3), 9L), "start id")
})

test_that("make_compiled_loo_parameters surfaces degeneracy as NA, not an error", {
  params <- make_compiled_loo_parameters(rep(0.5, 6), get_start("beta"))
  expect_true(anyNA(params))
  expect_equal(colnames(params), c("shape1", "shape2"))

  small <- make_compiled_loo_parameters(c(1, 2), get_start("gamma"))
  expect_true(all(is.na(small)))
})

test_that("compiled leave-one-out validation rejects non-finite parameters", {
  expect_false(compiled_loo_parameters_are_valid(
    matrix(Inf, nrow = 4L, ncol = 1L),
    2L
  ))
  expect_false(compiled_loo_parameters_are_valid(
    cbind(mean = 1:4, sd = 0),
    1L
  ))
})

test_that("exponential leave-one-out means avoid subtractive cancellation", {
  x <- c(rep(1e-300, 20), 1)
  parameters <- make_compiled_loo_parameters(x, get_start("exponential"))

  expect_true(all(is.finite(parameters)))
  expect_equal(
    unname(parameters[nrow(parameters), "rate"]),
    1e300,
    tolerance = 1e-12
  )
})

test_that("compiled leave-one-out fits handle numerical boundary cases", {
  exponential <- designated_loo_parameters_cpp(c(0, 0, 1), 2L)
  expect_true(is.na(exponential[3, 1]))

  almost_constant <- c(
    0x1.a2fe76a3f9475p-498,
    0x1.a2fe76a3f9475p-499,
    0x1.a2fe76a3f77aap-499,
    0x1.a2fe76a3f9473p-499,
    0x1.a2fe76a3f9473p-499,
    0x1.a2fe76a3fb14p-499,
    0x1.a2fe76a3f77aap-499,
    0x1.a2fe76a3f77aap-499,
    0x1.a2fe76a3f9473p-499
  )
  normal <- designated_loo_parameters_cpp(almost_constant, 1L)
  expect_identical(normal[1, 2], 0)

  weibull <- designated_loo_parameters_cpp(rep(1, 3), 5L)
  expect_true(all(is.na(weibull)))

  difficult_beta <- c(
    9.01349529317247e-16,
    8.34193070929293e-58,
    9.96483385752205e-90,
    1.3187543768896e-114,
    1.65819281622785e-189
  )
  beta <- designated_loo_parameters_cpp(difficult_beta, 6L)
  expect_true(all(is.na(beta)))
})

test_that("degenerate fits return the same NaN structure across designated starts", {
  ## n = 2 is degenerate for every designated start; all must return NaN, none
  ## may throw.
  for (start_id in 1:6) {
    m <- designated_loo_parameters_cpp(c(0.3, 0.6), start_id)
    expect_true(all(is.na(m)), info = paste("start_id", start_id))
  }
})

test_that("the compiled UCV term returns NA (not an error) on non-finite parameters", {
  x <- c(0.2, 0.4, 0.6, 0.8)
  h <- 0.1
  bad_params <- matrix(NA_real_, nrow = length(x), ncol = 2)
  expect_true(is.na(
    compiled_ucv_term2(x, h, get_kernel("beta"), get_start("beta"), bad_params)
  ))
})
