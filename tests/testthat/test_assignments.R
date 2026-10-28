# Fitted kdensity objects are read-only. There are no `$<-`, `[[<-` or
# update() methods; update() is stats::update.default(), which re-evaluates
# the stored call with the changed arguments, as for lm().

obj <- kdensity(cars$dist, kernel = "gamma")

test_that("a fresh kdensity object has the expected defaults", {
  expect_equal(obj$start_str, "uniform")
  expect_identical(getCall(obj), obj$call)
})

test_that("update() refits from the stored call without mutating its argument", {
  updated <- update(obj, start = "gamma")
  expect_equal(updated$start_str, "gamma")
  expect_equal(obj$start_str, "uniform")
  expect_equal(
    updated(c(10, 40)),
    kdensity(cars$dist, kernel = "gamma", start = "gamma")(c(10, 40))
  )
  expect_equal(
    update(obj, x = 1:100)(1:10),
    kdensity(1:100, kernel = "gamma")(1:10)
  )
})

test_that("update() keeps numeric and custom bandwidths", {
  set.seed(1)
  y <- rnorm(100)
  numeric_fit <- kdensity(y, bw = 0.5)
  updated <- update(numeric_fit, kernel = "epanechnikov")
  expect_equal(updated$kernel_str, "epanechnikov")
  expect_equal(updated$bw, 0.5)

  fixed_width <- function(x, kernel = NULL, start = NULL, support = NULL) 0.3
  custom_fit <- kdensity(y, bw = fixed_width)
  updated <- update(custom_fit, adjust = 2)
  expect_equal(updated$bw_str, "fixed_width")
  expect_equal(updated$h, 2 * custom_fit$h)
})

test_that("fitted objects cannot be modified in place", {
  fit <- obj
  expect_error(fit$start <- "gaussian")
  expect_error(fit[["adjust"]] <- 2)
  expect_equal(obj$start_str, "uniform")
})
