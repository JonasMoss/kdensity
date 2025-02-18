context("kernels")

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
