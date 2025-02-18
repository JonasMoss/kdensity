context("generics")

obj <- kdensity(mtcars$mpg, kernel = "gaussian", start = "gumbel")

expect_equal(obj$h, obj[["h"]])
expect_s3_class(logLik(obj), "logLik")

expect_equal(
  update(obj, x = 1:100)(1:10),
  kdensity(1:100, kernel = "gaussian", start = "gumbel")(1:10)
)
