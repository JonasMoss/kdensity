test_that("underflowing normal starts fail before bandwidth selection", {
  x <- c(rep(c(-1, 1), 5000), 1000)
  message <- paste0(
    "Parametric start 'normal'.*1 observation has.*x = 1000.*",
    "heavier-tailed start or check for outliers"
  )
  expect_error(kdensity(x, start = "normal"), message)
  expect_error(kdensity(x, start = "normal", bw = "nrd0"), message)
  never_called <- function(...) stop("bandwidth selector should not run")
  expect_error(kdensity(x, start = "normal", bw = never_called), message)
})

test_that("custom starts reject all non-positive and non-finite densities", {
  for (bad in c(0, -1, Inf, NaN, NA_real_)) {
    custom_start <- list(
      density = function(x) ifelse(x == 1, bad, 1),
      estimator = function(x) numeric(0),
      support = c(-Inf, Inf)
    )
    expect_error(
      kdensity(c(0, 1, 1), start = custom_start, bw = 1),
      "Parametric start 'custom_start'.*2 observations have.*x = 1"
    )
  }
})

test_that("numeric infinite bandwidth permits underflowing start densities", {
  x <- c(rep(c(-1, 1), 5000), 1000)
  fit <- kdensity(x, start = "normal", bw = Inf)
  y <- c(-1, 0, 1, 1000)
  expect_equal(
    fit(y),
    stats::dnorm(y, mean = fit$estimates[["mean"]], sd = fit$estimates[["sd"]])
  )
  expect_true(all(is.finite(fit(y))))
})

test_that("ordinary fits retain the weighted kernel estimate", {
  x <- c(-2, -1, 0, 1, 2)
  y <- c(-0.5, 0, 0.5)
  fit <- kdensity(x, start = "normal", bw = 0.5, normalized = FALSE)
  g <- function(z) stats::dnorm(
    z, mean = fit$estimates[["mean"]], sd = fit$estimates[["sd"]]
  )
  expected <- vapply(y, function(z) {
    g(z) * mean(stats::dnorm(z - x, sd = 0.5) / g(x))
  }, numeric(1))
  expect_equal(fit(y), expected)
  expect_equal(kdensity(x, bw = 0.5)(y), vapply(y, function(z) {
    mean(stats::dnorm(z - x, sd = 0.5))
  }, numeric(1)))
})
