#!/usr/bin/env Rscript

# Benchmark exact evaluation against the optional FFT grid backend.
#
# Rscript tools/grid-benchmark.R \
#   --reps=10 --sizes=100,1000,10000 --grid_n=1000 --fft_n=4096

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1]]))
repository <- dirname(dirname(script_path))

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("The benchmark requires the development package 'pkgload'.")
}
pkgload::load_all(repository, quiet = TRUE)

parse_options <- function(arguments) {
  options <- list(
    reps = 10L,
    sizes = c(100L, 1000L, 10000L),
    grid_n = 1000L,
    fft_n = 4096L,
    seed = 20260720L,
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
      options$grid_n < 2L || options$fft_n < 512L) {
    stop("Use reps >= 1, sizes >= 10, grid_n >= 2, and fft_n >= 512.")
  }
  options
}

median_elapsed <- function(fun, repetitions) {
  median(replicate(
    repetitions,
    unname(system.time(fun())[["elapsed"]])
  ))
}

make_case <- function(case, n) {
  switch(case,
    gaussian_uniform = kdensity::kdensity(
      stats::rnorm(n),
      kernel = "gaussian",
      start = "uniform",
      bw = 0.2
    ),
    gaussian_normal = kdensity::kdensity(
      stats::rnorm(n, mean = 0.4, sd = 1.3),
      kernel = "gaussian",
      start = "normal",
      bw = 0.2
    ),
    gcopula_uniform = kdensity::kdensity(
      stats::rbeta(n, 2, 5),
      kernel = "gcopula",
      start = "uniform",
      bw = 0.12
    ),
    gcopula_beta = kdensity::kdensity(
      stats::rbeta(n, 2, 5),
      kernel = "gcopula",
      start = "beta",
      bw = 0.12
    )
  )
}

benchmark_case <- function(case, n, options) {
  fit <- make_case(case, n)
  limits <- if (grepl("^gcopula", case)) {
    c(0.001, 0.999)
  } else {
    c(-3.5, 4.5)
  }
  fast <- kdensity::kdensity_grid(
    fit,
    from = limits[[1L]],
    to = limits[[2L]],
    n = options$grid_n,
    method = "fft",
    fft_n = options$fft_n
  )
  exact <- fit(fast$x)

  exact_seconds <- median_elapsed(
    function() fit(fast$x),
    options$reps
  )
  fft_seconds <- median_elapsed(
    function() kdensity::kdensity_grid(
      fit,
      from = limits[[1L]],
      to = limits[[2L]],
      n = options$grid_n,
      method = "fft",
      fft_n = options$fft_n
    ),
    options$reps
  )
  central <- exact > max(exact) * 1e-4

  data.frame(
    case = case,
    sample_size = n,
    grid_points = options$grid_n,
    fft_points = options$fft_n,
    exact_seconds = exact_seconds,
    fft_seconds = fft_seconds,
    speedup = exact_seconds / fft_seconds,
    max_absolute_error = max(abs(fast$y - exact)),
    root_mean_squared_error = sqrt(mean((fast$y - exact)^2)),
    max_central_relative_error = max(
      abs(fast$y[central] - exact[central]) / exact[central]
    )
  )
}

options <- parse_options(commandArgs(trailingOnly = TRUE))
set.seed(options$seed)
cases <- c(
  "gaussian_uniform",
  "gaussian_normal",
  "gcopula_uniform",
  "gcopula_beta"
)
results <- do.call(rbind, lapply(options$sizes, function(n) {
  do.call(rbind, lapply(cases, benchmark_case, n = n, options = options))
}))
row.names(results) <- NULL
print(results)

if (nzchar(options$output)) {
  directory <- dirname(options$output)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(results, options$output, row.names = FALSE)
  message("Wrote: ", options$output)
}
