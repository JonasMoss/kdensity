context("kdensity")

# Normalization errors.
x <- c(170, 172, 142, 199, 180, 184)
expect_error(kdensity(x, start = "normal"))
expect_error(kdensity(x, kernel = "epanechnikov"))
set.seed(1)
x <- rweibull(100, 2, 7)
expect_error(kdensity(x, kernel = "triangular", start = "weibull"))

expect_error(kdensity(c(precip, NA)), "x contains NAs and na.rm = FALSE.")
expect_equal(kdensity(c(precip, NA), na.rm = TRUE)(10:20), kdensity(precip)(10:20))
normal2 <- list(
  density = dnorm,
  estimator = function(data) {
    c(
      mean = mean(data),
      sd = sqrt(stats::var(data) * (length(data) - 1) / length(data))
    )
  },
  support = c(-Inf, Inf)
)

gaussian2 <- list(
  kernel  = function(y, x, h) dnorm((y - x) / h),
  sd      = 1,
  support = c(-Inf, Inf)
)

set.seed(313)
silly_width <- function(x, kernel = NULL, start = NULL, support = NULL) {
  1
}

expect_equal(
  kdensity(precip, start = "normal")(10),
  kdensity(precip, start = normal2)(10)
)
expect_equal(
  kdensity(precip, start = "normal", kernel = "gaussian")(10),
  kdensity(precip, start = normal2, kernel = gaussian2)(10)
)
expect_error(kdensity(precip, kernel = "beta"))
expect_error(kdensity(precip, support = c(0, 1)))
expect_error(kdensity(precip, bw = Inf))
expect_equal(
  kdensity(precip, bw = Inf, start = "normal")(10),
  dnorm(10,
    mean = mean(precip),
    sd = sqrt(stats::var(precip) * (length(precip) - 1) / length(precip))
  )
)
expect_equal(
  kdensity(precip, bw = 1)(10),
  kdensity(precip, bw = silly_width)(10)
)
expect_error(kdensity(precip)())
expect_error(kdensity(precip,
  start = "gumbel", kernel = "rectangular",
  bw = "RHE"
))
