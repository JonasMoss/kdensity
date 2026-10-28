test_that("sorted compact evaluators match direct compiled summation", {
  set.seed(20260721)
  kernel_ids <- 2:9

  for (offset in c(0, 1e8)) {
    x <- sort(offset + stats::rnorm(500, sd = 0.5))
    weights <- stats::runif(length(x))
    h <- 0.2
    y <- c(
      stats::setNames(
        offset + seq(-1.5, 1.5, length.out = 301),
        paste0("point", seq_len(301))
      ),
      missing = NA_real_,
      positive_infinity = Inf,
      negative_infinity = -Inf
    )
    cache <- compact_kernel_cache_cpp(x, weights, h)

    for (kernel_id in kernel_ids) {
      actual <- kernel_weighted_sum_sorted_cpp(
        y, x, weights, h, kernel_id, cache
      )
      expected <- kernel_weighted_sum_cpp(
        y, x, weights, h, kernel_id, numeric(), numeric()
      )

      expect_equal(
        actual,
        expected,
        tolerance = 1e-11,
        info = paste("kernel", kernel_id, "offset", offset)
      )
      expect_named(actual, names(y))
    }
  }
})

test_that("the sorted Laplace recurrence matches direct summation", {
  set.seed(20260721)
  x <- sort(c(
    stats::rnorm(500),
    rep(c(-0.25, 0.25), each = 3)
  ))
  weights <- stats::runif(length(x))
  y <- c(
    right = 1,
    duplicate1 = 0.25,
    left = -1,
    duplicate2 = 0.25,
    missing = NA_real_,
    positive_infinity = Inf
  )

  for (h in c(0.02, 0.2, 2)) {
    actual <- kernel_weighted_sum_sorted_cpp(
      y, x, weights, h, 10L, list()
    )
    expected <- kernel_weighted_sum_cpp(
      y, x, weights, h, 10L, numeric(), numeric()
    )

    expect_equal(actual, expected, tolerance = 1e-11)
    expect_named(actual, names(y))
  }
})

test_that("kernel evaluators switch after enough cumulative work", {
  set.seed(20260721)
  x <- stats::rnorm(1000)

  compact <- make_kernel_evaluator(
    x = x,
    h = 0.2,
    kernel = get_kernel("triweight"),
    parametric_start_data = rep(1, length(x))
  )
  compact(0)
  expect_false(environment(compact)$sorted_ready)
  compact(seq(-2, 2, length.out = 200))
  expect_true(environment(compact)$sorted_ready)

  laplace <- make_kernel_evaluator(
    x = x,
    h = 0.2,
    kernel = get_kernel("laplace"),
    parametric_start_data = rep(1, length(x))
  )
  laplace(seq(-1, 1, length.out = 9))
  expect_false(environment(laplace)$sorted_ready)
  laplace(c(-0.5, 0.5))
  expect_true(environment(laplace)$sorted_ready)
})

test_that("the dispatched sorted path retains parametric-start weights", {
  set.seed(20260721)
  x <- stats::rnorm(1000)
  start_data <- stats::dnorm(x, mean = 0.2, sd = 1.3)
  y <- seq(-2, 2, length.out = 251)

  for (kernel_name in c(
    "epanechnikov", "rectangular", "triangular", "biweight",
    "cosine", "optcosine", "triweight", "tricube", "laplace"
  )) {
    kernel <- get_kernel(kernel_name)
    evaluator <- make_kernel_evaluator(
      x = x,
      h = 0.2,
      kernel = kernel,
      parametric_start_data = start_data
    )
    actual <- evaluator(y)
    expected <- vapply(y, function(point) {
      mean(kernel$kernel(point, x, 0.2) / start_data) / 0.2
    }, numeric(1))

    expect_equal(
      actual,
      expected,
      tolerance = 1e-10,
      info = kernel_name
    )
  }
})

test_that("sorted-kernel helpers validate their inputs", {
  expect_error(
    kernel_weighted_sum_sorted_cpp(
      0, numeric(), numeric(), 0.2, 2L, list()
    ),
    "equal positive length"
  )
  expect_error(
    compact_kernel_cache_cpp(c(1, 0), c(1, 1), 0.2),
    "must be sorted"
  )
  expect_error(
    compact_kernel_cache_cpp(c(0, 1), 1, 0.2),
    "equal positive length"
  )
  expect_error(
    compact_kernel_cache_cpp(c(0, 1), c(1, 1), 0),
    "finite and positive"
  )
  expect_error(
    compact_kernel_cache_cpp(c(0, Inf), c(1, 1), 0.2),
    "must be finite"
  )
  expect_error(
    kernel_weighted_sum_sorted_cpp(
      0, c(0, 1), c(1, 1), 0.2, 1L, list()
    ),
    "kernel ids 2 through 10"
  )
  expect_error(
    kernel_weighted_sum_sorted_cpp(
      0, c(0, 1), c(1, 1), 0, 2L, list()
    ),
    "finite and positive"
  )
  expect_error(
    kernel_weighted_sum_sorted_cpp(
      0, c(0, Inf), c(1, 1), 0.2, 2L, list()
    ),
    "must be finite"
  )
  expect_error(
    kernel_weighted_sum_sorted_cpp(
      0, c(1, 0), c(1, 1), 0.2, 2L, list()
    ),
    "must be sorted"
  )
  expect_error(
    kernel_weighted_sum_sorted_cpp(
      0, c(0, 1), c(1, 1), 0.2, 2L, list()
    ),
    "cache is incomplete"
  )
  expect_error(
    kernel_weighted_sum_sorted_cpp(
      0, c(0, 1), c(1, 1), 0.2, 2L,
      list(
        base = 2L,
        centers = numeric(3),
        summaries = matrix(0, nrow = 4, ncol = 14)
      )
    ),
    "incompatible dimensions"
  )
})

test_that("uncached compact nodes recurse to their leaves", {
  cache <- list(
    base = 2L,
    centers = numeric(4),
    summaries = matrix(NA_real_, nrow = 4, ncol = 14)
  )
  expect_equal(
    kernel_weighted_sum_sorted_cpp(
      0.5, c(0, 1), c(1, 1), 1, 2L, cache
    ),
    0
  )
})

test_that("the sorted rectangular kernel excludes boundary ties like the direct path", {
  # Regression: the sorted path selected its window with subtraction bounds
  # (point +/- h), but the direct path and stats::density include x iff
  # abs((x - point) / h) <= 1. These differ by one ULP at the boundary; the
  # rectangular kernel jumps to 0.5 there, so an excluded boundary tie used to
  # contribute in full (0.4286 instead of 0 for the case below). Continuous
  # random data never lands on the boundary, so the other tests miss this.
  x <- sort(c(rep(1.0, 20), rep(2.0, 60), rep(3.5, 20)))
  weights <- rep(1 / length(x), length(x))
  h <- 0.7
  y <- c(boundary = 2.7, seq(0, 5, length.out = 101))

  cache <- compact_kernel_cache_cpp(x, weights, h)
  actual <- kernel_weighted_sum_sorted_cpp(y, x, weights, h, 3L, cache)
  expected <- kernel_weighted_sum_cpp(y, x, weights, h, 3L, numeric(), numeric())

  expect_equal(actual, expected, tolerance = 1e-12)
  # (2.0 - 2.7) / 0.7 == -1.0000000000000002 -> excluded -> exact value is 0.
  expect_equal(unname(actual[["boundary"]]), 0)
})

test_that("every sorted compact kernel matches direct summation on boundary ties", {
  set.seed(20260721)
  scenarios <- list(
    list(x = c(rep(1, 20), rep(2, 60), rep(3.5, 20)), h = 0.7),
    list(x = rep(seq(0, 3, by = 0.5), each = 20), h = 0.5),
    list(x = rep(c(0.1, 0.4, 0.4, 0.9, 1.6), each = 25), h = 0.3)
  )
  for (s in scenarios) {
    x <- sort(as.double(s$x))
    weights <- stats::runif(length(x))
    ## Evaluation points include exact window edges (x +/- h) and interior ties.
    y <- sort(unique(c(x, x + s$h, x - s$h, x + s$h / 2)))
    cache <- compact_kernel_cache_cpp(x, weights, s$h)
    for (kernel_id in 2:10) {
      cache_used <- if (kernel_id == 10L) list() else cache
      actual <- kernel_weighted_sum_sorted_cpp(
        y, x, weights, s$h, kernel_id, cache_used
      )
      expected <- kernel_weighted_sum_cpp(
        y, x, weights, s$h, kernel_id, numeric(), numeric()
      )
      expect_equal(actual, expected, tolerance = 1e-12,
        info = paste("kernel", kernel_id, "h", s$h)
      )
    }
  }
})
