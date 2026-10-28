t_without_kernel <- list(
  support = c(-Inf, Inf)
)

t_without_support <- list(
  kernel = dt
)

t <- list(
  kernel = dt,
  support = c(-Inf, Inf)
)

expect_equal(get_kernel("normal")$sd, 1)
expect_silent(add_kernel("t", t))
expect_error(add_kernel("t", t_without_kernel))
expect_error(add_kernel("t", t_without_support))
expect_error(get_kernel(5))
expect_error(get_kernel("no_kernel"))

test_that("the gamma kernels follow Chen (2000)", {
  # Written from the paper rather than from the package code: the first kernel
  # is K_{y/b + 1, b} and the modified kernel is K_{rho_b(y), b}, where
  # rho_b(y) = y / b for y >= 2b and (y / (2b))^2 + 1 for 0 <= y < 2b.
  chen_rho <- function(y, b) ifelse(y >= 2 * b, y / b, (y / (2 * b))^2 + 1)
  chen_modified <- function(y, t, b) dgamma(t, shape = chen_rho(y, b), scale = b)
  chen_first <- function(y, t, b) dgamma(t, shape = y / b + 1, scale = b)

  set.seed(20261028)
  b <- 0.15
  t <- c(0, rgamma(40, 2, 1.5))
  y <- c(0, 0.5 * b, b, 2 * b - 1e-12, 2 * b, 2.5 * b, 1, 3)
  modified <- get_kernel("gamma")
  first <- get_kernel("gamma_biased")
  modified_evaluator <- make_kernel_evaluator(t, b, modified, rep(1, length(t)))
  first_evaluator <- make_kernel_evaluator(t, b, first, rep(1, length(t)))

  for (point in y) {
    expect_equal(modified$kernel(point, t, b) / b, chen_modified(point, t, b),
      tolerance = 1e-12, info = paste("R modified, y =", point)
    )
    expect_equal(first$kernel(point, t, b) / b, chen_first(point, t, b),
      tolerance = 1e-12, info = paste("R first, y =", point)
    )
  }
  expect_equal(modified_evaluator(y),
    vapply(y, function(point) mean(chen_modified(point, t, b)), numeric(1)),
    tolerance = 1e-12
  )
  expect_equal(first_evaluator(y),
    vapply(y, function(point) mean(chen_first(point, t, b)), numeric(1)),
    tolerance = 1e-12
  )

  # In the interior the modified kernel has mean y and the first kernel y + b,
  # so the two kernels must differ there.
  kernel_mean <- function(kernel, point) {
    integrate(function(s) s * kernel$kernel(point, s, b) / b, 0, Inf)$value
  }
  expect_equal(kernel_mean(modified, 1), 1, tolerance = 1e-8)
  expect_equal(kernel_mean(first, 1), 1 + b, tolerance = 1e-8)
  # The modified shape is continuous at the cut point 2b.
  expect_equal(chen_rho(2 * b - 1e-12, b), chen_rho(2 * b, b), tolerance = 1e-9)
})
