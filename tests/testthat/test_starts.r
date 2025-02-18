context("parametric starts")

t_without_estim <- list(
  density = dt,
  support = c(-Inf, Inf)
)

t_without_support <- list(
  density = dt,
  estimator = function(x) NULL
)

t_without_density <- list(
  estimator = function(x) NULL,
  support = c(-Inf, Inf)
)

t <- list(
  density = dt,
  estimator = function(x) NULL,
  support = c(-Inf, Inf)
)

expect_equal(get_start("normal")$support, c(-Inf, Inf))
expect_silent(add_start("t", t))
expect_error(add_start("t", t_without_estim))
expect_error(add_start("t", t_without_support))
expect_error(add_start("t", t_without_density))
expect_error(get_start(5))
expect_error(get_start("no_start"))
