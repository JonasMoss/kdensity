context("assignments and updates")

cars$dist %>%
  kdensity(kernel = "gamma") ->
obj

expect_equal(obj$start_str, "uniform")

update(obj, start = "gamma")
expect_equal(obj$start_str, "gamma")

obj[["start"]] <- "gaussian"
expect_equal(obj$start_str, "gaussian")

obj$start <- "gumbel"
expect_equal(obj$start_str, "gumbel")
