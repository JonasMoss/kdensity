#!/usr/bin/env Rscript

# Benchmark the exact sorted-data evaluator against direct O(n * m)
# summation. The "first_sorted" timing includes sorting and cache construction;
# "cached_sorted" measures later evaluations of the same fitted object.
#
# Rscript tools/exact-kernel-benchmark.R \
#   --reps=5 --sizes=100,1000,10000 --points=1000

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1]]))
repository <- dirname(dirname(script_path))

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("The benchmark requires the development package 'pkgload'.")
}
pkgload::load_all(repository, quiet = TRUE)

parse_options <- function(arguments) {
  options <- list(
    reps = 5L,
    sizes = c(100L, 1000L, 10000L),
    points = 1000L,
    seed = 20260721L,
    output = ""
  )

  for (argument in arguments) {
    pieces <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1]]
    if (length(pieces) != 2L || !pieces[[1]] %in% names(options)) {
      stop("Unknown option: ", argument)
    }
    name <- pieces[[1]]
    value <- pieces[[2]]
    if (name == "sizes") {
      options[[name]] <- as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
    } else if (name == "output") {
      options[[name]] <- value
    } else {
      options[[name]] <- as.integer(value)
    }
  }

  if (options$reps < 1L || any(options$sizes < 10L) ||
      options$points < 2L) {
    stop("Use reps >= 1, sizes >= 10, and points >= 2.")
  }
  options
}

median_seconds <- function(fun, repetitions, batch) {
  median(replicate(repetitions, {
    start <- proc.time()[["elapsed"]]
    for (iteration in seq_len(batch)) {
      fun()
    }
    (proc.time()[["elapsed"]] - start) / batch
  }))
}

benchmark_case <- function(kernel_name, n, options) {
  kernel <- kdensity:::get_kernel(kernel_name)
  kernel_id <- kernel$compiled_id
  h <- 0.2
  x <- stats::rnorm(n)
  y <- seq(-2.5, 2.5, length.out = options$points)
  weights <- rep(1 / (n * h), n)
  empty_cache <- numeric()
  batch <- max(
    1L,
    as.integer(ceiling(1e7 / (as.double(n) * options$points)))
  )

  direct <- function() {
    kdensity:::kernel_weighted_sum_cpp(
      y, x, weights, h, kernel_id, empty_cache, empty_cache
    )
  }
  first_sorted <- function() {
    evaluator <- kdensity:::make_kernel_evaluator(
      x, h, kernel, rep(1, n)
    )
    evaluator(y)
  }
  cached_evaluator <- kdensity:::make_kernel_evaluator(
    x, h, kernel, rep(1, n)
  )
  cached_evaluator(y)
  cached_sorted <- function() cached_evaluator(y)

  expected <- direct()
  actual <- cached_sorted()
  direct_seconds <- median_seconds(direct, options$reps, batch)
  first_seconds <- median_seconds(
    first_sorted, options$reps, batch
  )
  cached_seconds <- median_seconds(
    cached_sorted, options$reps, batch
  )

  data.frame(
    kernel = kernel_name,
    sample_size = n,
    evaluation_points = options$points,
    direct_seconds = direct_seconds,
    first_sorted_seconds = first_seconds,
    cached_sorted_seconds = cached_seconds,
    first_speedup = direct_seconds / first_seconds,
    cached_speedup = direct_seconds / cached_seconds,
    max_absolute_error = max(abs(actual - expected)),
    root_mean_squared_error = sqrt(mean((actual - expected)^2))
  )
}

options <- parse_options(commandArgs(trailingOnly = TRUE))
set.seed(options$seed)
kernels <- c(
  "epanechnikov", "rectangular", "triangular", "biweight",
  "cosine", "optcosine", "triweight", "tricube", "laplace"
)
results <- do.call(rbind, lapply(options$sizes, function(n) {
  do.call(rbind, lapply(
    kernels,
    benchmark_case,
    n = n,
    options = options
  ))
}))
row.names(results) <- NULL
print(results)

if (nzchar(options$output)) {
  directory <- dirname(options$output)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(results, options$output, row.names = FALSE)
  message("Wrote: ", options$output)
}
