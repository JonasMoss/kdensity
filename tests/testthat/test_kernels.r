context("kernels")

t_without_kernel = list(
  support = c(-Inf, Inf)
)

t_without_support = list(
  kernel = dt
)

expect_error(add_kernel("t", t_without_kernel))
expect_error(add_kernel("t", t_without_support))
expect_error(get_kernel(5))
expect_error(get_kernel("no_kernel"))
