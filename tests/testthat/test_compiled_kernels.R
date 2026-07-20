test_that("compiled built-in kernels match their R definitions", {
  set.seed(42)

  cases <- list(
    gaussian = list(x = rnorm(31), y = seq(-2, 2, length.out = 17), h = 0.3),
    epanechnikov = list(x = rnorm(31), y = seq(-2, 2, length.out = 17), h = 0.3),
    rectangular = list(x = rnorm(31), y = seq(-2, 2, length.out = 17), h = 0.3),
    triangular = list(x = rnorm(31), y = seq(-2, 2, length.out = 17), h = 0.3),
    biweight = list(x = rnorm(31), y = seq(-2, 2, length.out = 17), h = 0.3),
    cosine = list(x = rnorm(31), y = seq(-2, 2, length.out = 17), h = 0.3),
    optcosine = list(x = rnorm(31), y = seq(-2, 2, length.out = 17), h = 0.3),
    triweight = list(x = rnorm(31), y = seq(-2, 2, length.out = 17), h = 0.3),
    tricube = list(x = rnorm(31), y = seq(-2, 2, length.out = 17), h = 0.3),
    laplace = list(x = rnorm(31), y = seq(-2, 2, length.out = 17), h = 0.3),
    gcopula = list(x = rbeta(31, 2, 5), y = seq(0.01, 0.99, length.out = 17), h = 0.15),
    gamma = list(x = c(0, rgamma(30, 2, 2)), y = seq(0, 3, length.out = 17), h = 0.08),
    gamma_biased = list(x = c(0, rgamma(30, 2, 2)), y = seq(0, 3, length.out = 17), h = 0.08),
    beta = list(x = c(0, rbeta(29, 2, 5), 1), y = seq(0, 1, length.out = 17), h = 0.02),
    beta_biased = list(x = c(0, rbeta(29, 2, 5), 1), y = seq(0, 1, length.out = 17), h = 0.02)
  )

  for (kernel_name in names(cases)) {
    case <- cases[[kernel_name]]
    kernel <- get_kernel(kernel_name)
    parametric_start_data <- seq(0.75, 1.25, length.out = length(case$x))
    evaluator <- make_kernel_evaluator(
      x = case$x,
      h = case$h,
      kernel = kernel,
      parametric_start_data = parametric_start_data
    )

    expected <- vapply(case$y, function(point) {
      mean(
        kernel$kernel(point, case$x, case$h) /
          parametric_start_data
      ) / case$h
    }, numeric(1))

    expect_equal(
      evaluator(case$y),
      expected,
      tolerance = 1e-11,
      info = kernel_name
    )
  }
})

test_that("compiled evaluators preserve names and missing values", {
  x <- c(-1, 0, 1)
  y <- c(left = -0.5, missing = NA_real_, right = 0.5)
  evaluator <- make_kernel_evaluator(
    x = x,
    h = 0.4,
    kernel = get_kernel("gaussian"),
    parametric_start_data = rep(1, length(x))
  )

  result <- evaluator(y)
  expect_named(result, names(y))
  expect_true(is.na(result[["missing"]]))
})

test_that("custom kernels continue to use the R fallback", {
  calls <- 0L
  custom <- list(
    kernel = function(y, x, h) {
      calls <<- calls + 1L
      dnorm((y - x) / h)
    },
    support = c(-Inf, Inf),
    sd = 1
  )
  x <- c(-1, 0, 1)
  h <- 0.5
  start_data <- c(0.8, 1, 1.2)
  evaluator <- make_kernel_evaluator(
    x = x,
    h = h,
    kernel = custom,
    parametric_start_data = start_data
  )

  y <- c(-0.5, 0.5)
  weights <- 1 / (length(x) * h * start_data)
  expected <- vapply(y, function(point) {
    sum(dnorm((point - x) / h) * weights)
  }, numeric(1))

  ## The R fallback must reproduce the weighted kernel sum exactly, one kernel
  ## call per evaluation point.
  expect_equal(evaluator(y), expected)
  expect_equal(calls, length(y))
})

test_that("known built-in normalizations avoid numerical integration", {
  uniform <- get_start("uniform")
  normal <- get_start("normal")
  expect_equal(
    known_normalization(
      uniform, get_kernel("gaussian"), c(-Inf, Inf), 0.5,
      x = c(-1, 1), estimates = numeric()
    ),
    1
  )
  expect_equal(
    known_normalization(
      uniform, get_kernel("gcopula"), c(0, 1), 0.2,
      x = c(0.25, 0.75), estimates = numeric()
    ),
    sqrt(1 - (1 - 0.2^2)^2) / 0.2
  )
  expect_null(
    known_normalization(
      uniform, get_kernel("beta"), c(0, 1), 0.02,
      x = c(0.25, 0.75), estimates = numeric()
    )
  )
  expect_null(
    known_normalization(
      normal, get_kernel("beta"), c(0, 1), 0.02,
      x = c(0.25, 0.75), estimates = c(mean = 0.5, sd = 0.25)
    )
  )
})

test_that("normal starts with Gaussian kernels have a closed normalization", {
  x <- c(-1.5, -0.25, 0.5, 2)
  estimates <- c(mean = mean(x), sd = sqrt(mean((x - mean(x))^2)))
  h <- 0.35
  start_density <- dnorm(x, estimates[["mean"]], estimates[["sd"]])

  actual <- known_normalization(
    get_start("normal"),
    get_kernel("gaussian"),
    c(-Inf, Inf),
    h,
    x,
    estimates
  )
  expected <- integrate(function(y) {
    vapply(y, function(point) {
      mean(
        dnorm((point - x) / h) / (h * start_density)
      ) * dnorm(point, estimates[["mean"]], estimates[["sd"]])
    }, numeric(1))
  }, -Inf, Inf)$value

  expect_equal(actual, expected, tolerance = 1e-10)
})

test_that("embedded leave-one-out MLEs match designated starts", {
  set.seed(20260719)
  samples <- list(
    norm = rnorm(30, 0.4, 1.2),
    exp = rexp(30, 1.3),
    lnorm = rlnorm(30, 0.2, 0.7),
    gamma = rgamma(30, 2.8, 1.4),
    weibull = rweibull(30, 2.2, 1.5),
    beta = rbeta(30, 2.5, 4.1)
  )

  for (start_name in names(samples)) {
    x <- samples[[start_name]]
    start <- get_start(start_name)
    parameter_names <- compiled_start_parameter_names(start$compiled_id)
    actual <- make_compiled_loo_parameters(x, start)
    expected <- matrix(
      unlist(lapply(seq_along(x), function(i) {
        unname(start$estimator(x[-i])[parameter_names])
      })),
      nrow = length(x),
      byrow = TRUE
    )

    expect_equal(
      unname(actual),
      expected,
      tolerance = 1e-10,
      info = start_name
    )
  }
})

test_that("compiled non-uniform UCV terms match the R criterion", {
  set.seed(20260719)
  samples <- list(
    norm = rnorm(30, 0.4, 1.2),
    exp = rexp(30, 1.3),
    gamma = rgamma(30, 2.8, 1.4),
    weibull = rweibull(30, 2.2, 1.5),
    beta = rbeta(30, 2.5, 4.1)
  )
  cases <- list(
    c("gaussian", "norm"),
    c("gamma", "exp"),
    c("gamma", "gamma"),
    c("gamma", "weibull"),
    c("beta", "beta"),
    c("gcopula", "beta")
  )

  for (case in cases) {
    kernel_name <- case[[1]]
    start_name <- case[[2]]
    x <- samples[[start_name]]
    h <- if (kernel_name %in% c("beta", "gcopula")) 0.08 else 0.25
    kernel <- get_kernel(kernel_name)
    start <- get_start(start_name)
    parameters <- make_compiled_loo_parameters(x, start)
    x_name <- names(formals(start$density))[[1]]

    dstart <- function(data, estimates) {
      arguments <- c(list(data), as.list(estimates))
      names(arguments)[[1]] <- x_name
      do.call(start$density, arguments)
    }
    expected <- 2 * mean(vapply(seq_along(x), function(i) {
      keep <- seq_along(x) != i
      # Eval-first: shape from the eval point x[i], evaluated at x[keep].
      mean(
        kernel$kernel(x[i], x[keep], h) / h /
          dstart(x[keep], parameters[i, ])
      ) * dstart(x[i], parameters[i, ])
    }, numeric(1)))

    expect_equal(
      compiled_ucv_term2(x, h, kernel, start, parameters),
      expected,
      tolerance = 1e-11,
      info = paste(kernel_name, start_name)
    )
  }
})

test_that("the bias-corrected beta kernel preserves vector order", {
  kernel <- get_kernel("beta")$kernel
  y <- c(upper = 0.9, lower = 0.1, interior = 0.5)

  expect_equal(
    unname(kernel(y, x = 0.4, h = 0.08)),
    unname(vapply(y, kernel, numeric(1), x = 0.4, h = 0.08))
  )
})
