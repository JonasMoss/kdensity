# High-signal, deterministic invariants across many kernel/start combinations.

symmetric_kernels <- c(
  "gaussian", "epanechnikov", "rectangular", "triangular", "biweight",
  "triweight", "tricube", "cosine", "optcosine", "laplace"
)

test_that("built-in densities are finite and non-negative", {
  set.seed(20260720)
  cases <- list(
    list(x = rnorm(150), kernel = "gaussian", start = "uniform"),
    list(x = rnorm(150), kernel = "epanechnikov", start = "normal"),
    list(x = rnorm(150), kernel = "triangular", start = "gaussian"),
    list(x = abs(rnorm(150)) + 0.1, kernel = "gamma", start = "uniform"),
    list(x = rgamma(150, 3, 2), kernel = "gamma", start = "gamma"),
    list(x = rgamma(150, 3, 2), kernel = "gamma_biased", start = "gamma"),
    list(x = rbeta(150, 2, 5), kernel = "beta", start = "uniform"),
    list(x = rbeta(150, 2, 5), kernel = "beta_biased", start = "uniform"),
    list(x = rbeta(150, 2, 5), kernel = "gcopula", start = "uniform")
  )
  for (case in cases) {
    fit <- kdensity(case$x, kernel = case$kernel, start = case$start)
    grid <- seq(fit$range[1], fit$range[2], length.out = 300)
    grid <- grid[grid > fit$support[1] & grid < fit$support[2]]
    vals <- fit(grid)
    label <- paste(case$kernel, case$start)
    expect_true(all(is.finite(vals)), info = label)
    expect_true(all(vals >= -1e-9), info = label)
  }
})

test_that("closed-form-normalized densities integrate to one", {
  set.seed(20260720)
  x <- rnorm(200)
  for (kernel_name in c("gaussian", "epanechnikov", "triangular", "cosine")) {
    fit <- kdensity(x, kernel = kernel_name, start = "uniform")
    expect_equal(integrate(fit, -Inf, Inf)$value, 1,
      tolerance = 1e-3, info = kernel_name
    )
  }
  fit_n <- kdensity(x, kernel = "gaussian", start = "normal")
  expect_equal(integrate(fit_n, -Inf, Inf)$value, 1, tolerance = 1e-3)
})

test_that("numerically-normalized asymmetric densities integrate to about one", {
  set.seed(20260720)
  fg <- kdensity(rgamma(200, 3, 2), kernel = "gamma", start = "uniform")
  expect_equal(integrate(fg, 0, Inf)$value, 1, tolerance = 1e-2)

  fb <- kdensity(rbeta(200, 2, 5), kernel = "beta", start = "uniform")
  expect_equal(integrate(fb, 0, 1)$value, 1, tolerance = 1e-2)
})

test_that("symmetric-kernel estimates are reflection symmetric", {
  set.seed(20260720)
  x <- rnorm(120)
  grid <- seq(-3, 3, length.out = 61)
  for (kernel_name in symmetric_kernels) {
    f1 <- kdensity(x, kernel = kernel_name, start = "uniform", bw = 0.35)
    f2 <- kdensity(-x, kernel = kernel_name, start = "uniform", bw = 0.35)
    expect_equal(f1(grid), f2(-grid), tolerance = 1e-10, info = kernel_name)
  }
})

test_that("the compiled direct evaluator matches the R kernel definition", {
  set.seed(20260720)
  x <- rnorm(80)
  h <- 0.3
  start_data <- dnorm(x, 0.1, 1.1)
  weights <- 1 / (length(x) * h * start_data)
  y <- c(seq(-3, 3, length.out = 51), NA_real_)
  for (kernel_name in symmetric_kernels) {
    kernel <- get_kernel(kernel_name)
    evaluator <- make_kernel_evaluator(x, h, kernel, start_data)
    expected <- vapply(y, function(point) {
      sum(kernel$kernel(point, x, h) * weights)
    }, numeric(1))
    expect_equal(evaluator(y), expected, tolerance = 1e-9, info = kernel_name)
  }
})

test_that("non-uniform asymmetric UCV matches the eval-first leave-one-out term", {
  ds_of <- function(start) {
    nm <- names(formals(start$density))[[1]]
    function(data, est) {
      arguments <- c(list(data), as.list(est))
      names(arguments)[[1]] <- nm
      do.call(start$density, arguments)
    }
  }
  gold_term <- function(smp, kernel, start, h) {
    params <- make_compiled_loo_parameters(smp, start)
    ds <- ds_of(start)
    2 * mean(vapply(seq_along(smp), function(i) {
      ev <- make_kernel_evaluator(smp[-i], h, kernel, ds(smp[-i], params[i, ]))
      ev(smp[i]) * ds(smp[i], params[i, ])
    }, numeric(1)))
  }

  set.seed(7)
  cases <- list(
    list(k = "gamma", s = "gamma", x = rgamma(40, 3, 2), h = 0.25),
    list(k = "beta", s = "beta", x = rbeta(40, 2, 4), h = 0.08)
  )
  for (case in cases) {
    kernel <- get_kernel(case$k)
    start <- get_start(case$s)
    params <- make_compiled_loo_parameters(case$x, start)
    compiled <- compiled_ucv_term2(case$x, case$h, kernel, start, params)
    expect_equal(compiled, gold_term(case$x, kernel, start, case$h),
      tolerance = 1e-10, info = paste(case$k, case$s)
    )
  }
})

test_that("coef, logLik, AIC and BIC agree with the fitted parametric start", {
  set.seed(20260720)
  x <- rgamma(120, 3, 2)
  fit <- kdensity(x, start = "gamma", kernel = "gaussian")

  est <- coef(fit)
  expect_equal(as.numeric(est), as.numeric(get_start("gamma")$estimator(x)),
    tolerance = 1e-8
  )

  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  ## logLik is the log-likelihood of the fitted parametric start -- sum of the
  ## LOG densities, not the sum of densities.
  expect_equal(as.numeric(ll),
    sum(dgamma(x, shape = est[["shape"]], rate = est[["rate"]], log = TRUE)),
    tolerance = 1e-8
  )
  expect_equal(attr(ll, "df"), 2L)
  expect_equal(attr(ll, "nobs"), length(x))

  expect_equal(AIC(fit), -2 * as.numeric(ll) + 2 * 2, tolerance = 1e-8)
  expect_equal(BIC(fit), -2 * as.numeric(ll) + log(length(x)) * 2,
    tolerance = 1e-8
  )
})

test_that("logLik is unavailable for uniform and constant starts", {
  expect_error(logLik(kdensity(rnorm(50), start = "uniform")), "non-uniform")
  expect_true(is.na(kdensity(rnorm(50), start = "uniform")$logLik))
})
