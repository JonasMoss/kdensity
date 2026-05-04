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

## Other settable parameters round-trip through update().
obj[["kernel"]] <- "gaussian"
expect_equal(obj$kernel_str, "gaussian")

obj[["adjust"]] <- 2
expect_equal(obj$adjust, 2)

obj[["bw"]] <- "nrd0"
expect_equal(obj$bw_str, "nrd0")

obj[["normalized"]] <- FALSE
expect_false(obj$normalized)
