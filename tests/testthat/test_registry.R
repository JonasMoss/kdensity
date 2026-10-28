test_that("custom objects cannot replace any built-in kernel or start", {
  x <- c(-1, -0.5, 0, 0.5, 1)
  kernel_before <- lapply(builtin_kernel_names, get_kernel)
  start_before <- lapply(builtin_start_names, get_start)
  values_before <- kdensity(x, kernel = "gaussian", start = "normal", bw = 0.5)(x)

  gaussian <- list(
    kernel = function(y, x, h) 2 * stats::dnorm((y - x) / h),
    sd = 1, support = c(-Inf, Inf)
  )
  normal <- list(
    density = stats::dnorm,
    estimator = function(data) c(mean = mean(data), sd = stats::sd(data)),
    support = c(-Inf, Inf)
  )
  for (name in builtin_kernel_names) {
    expect_error(add_kernel(name, gaussian), "built-in kernel.*Rename your object")
  }
  for (name in builtin_start_names) {
    expect_error(add_start(name, normal), "built-in parametric start.*Rename your object")
  }
  expect_error(kdensity(x, kernel = gaussian), "my_gaussian")
  expect_error(kdensity(x, start = normal), "my_normal")
  expect_identical(lapply(builtin_kernel_names, get_kernel), kernel_before)
  expect_identical(lapply(builtin_start_names, get_start), start_before)
  expect_identical(get_kernel("gaussian")$compiled_id, 1L)
  expect_identical(get_start("normal")$compiled_id, 1L)
  expect_identical(
    kdensity(x, kernel = "gaussian", start = "normal", bw = 0.5)(x),
    values_before
  )
})

test_that("custom kernels and starts work and can be re-registered", {
  on.exit({
    rm(list = "registry_kernel", envir = kernel_environment)
    rm(list = "registry_start", envir = starts_environment)
  })
  x <- c(-1, -0.5, 0, 0.5, 1)
  registry_kernel <- list(
    kernel = function(y, x, h) stats::dnorm((y - x) / h),
    sd = 1, support = c(-Inf, Inf)
  )
  registry_start <- list(
    density = stats::dnorm,
    estimator = function(data) c(mean = mean(data), sd = stats::sd(data)),
    support = c(-Inf, Inf)
  )
  for (iteration in 1:2) {
    custom <- kdensity(x, kernel = registry_kernel, start = registry_start, bw = 0.5)
    parameters <- as.list(registry_start$estimator(x))
    density <- function(y) do.call(registry_start$density, c(list(y), parameters))
    h <- 0.5 * registry_kernel$sd
    raw <- function(y) vapply(y, function(point) {
      density(point) * mean(registry_kernel$kernel(point, x, h) / density(x)) / h
    }, numeric(1))
    normalizer <- stats::integrate(raw, -Inf, Inf)$value
    expect_equal(custom(x), raw(x) / normalizer, tolerance = 1e-6)
    expect_equal(kdensity(x, kernel = "registry_kernel", start = "registry_start", bw = 0.5)(x),
      custom(x))
    expect_identical(get_kernel("registry_kernel"), registry_kernel)
    expect_identical(get_start("registry_start"), registry_start)
    registry_kernel$sd <- 2
    registry_start$estimator <- function(data) c(mean = 0, sd = 1)
    x <- c(-1, 0, 1)
  }
})
