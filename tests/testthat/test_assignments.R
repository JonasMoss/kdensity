# Assignment, update, and replacement behaviour for kdensity objects.
#
# `update()` follows R's copy-on-modify semantics: it returns a new object and
# never mutates the caller's binding. The `$<-` and `[[<-` methods pass that
# returned object back through R's replacement protocol.

obj <- kdensity(cars$dist, kernel = "gamma")

test_that("a fresh kdensity object has the expected defaults", {
  expect_equal(obj$start_str, "uniform")
})

test_that("update() returns a new object without mutating its argument", {
  updated <- update(obj, start = "gamma")
  expect_equal(updated$start_str, "gamma")
  ## Regression test: the original binding must be untouched. The previous
  ## implementation reassigned it in place via assign(..., parent.frame()).
  expect_equal(obj$start_str, "uniform")
})

test_that("update() rebuilds the estimate from the supplied arguments", {
  expect_equal(
    update(obj, x = 1:100)(1:10),
    kdensity(1:100, kernel = "gamma")(1:10)
  )
})

test_that("[[<- rebuilds the object and round-trips settable fields", {
  fit <- obj
  fit[["start"]] <- "gaussian"
  expect_equal(fit$start_str, "gaussian")

  fit[["kernel"]] <- "gaussian"
  expect_equal(fit$kernel_str, "gaussian")

  fit[["adjust"]] <- 2
  expect_equal(fit$adjust, 2)

  fit[["bw"]] <- "nrd0"
  expect_equal(fit$bw_str, "nrd0")

  fit[["normalized"]] <- FALSE
  expect_false(fit$normalized)

  ## Replacing fields on `fit` must not reach back and mutate `obj`.
  expect_equal(obj$start_str, "uniform")
})

test_that("$<- rebuilds the object like [[<-", {
  fit <- obj
  fit$start <- "gumbel"
  expect_equal(fit$start_str, "gumbel")
})

test_that("replacement functions compose on nested structures", {
  ## The old parent.frame() assignment made this depend on the object being
  ## bound to a bare symbol; it now works through the replacement protocol.
  fits <- list(obj)
  fits[[1]]$start <- "gaussian"
  expect_equal(fits[[1]]$start_str, "gaussian")
})
