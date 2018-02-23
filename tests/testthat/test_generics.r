context("generics")

obj = kdensity(mtcars$mpg, kernel = "gaussian", start = "gumbel")

expect_equal(obj$h, obj[["h"]])
expect_s3_class(logLik(obj), "logLik")
