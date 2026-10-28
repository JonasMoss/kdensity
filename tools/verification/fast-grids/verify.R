#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1]]))
verification_dir <- dirname(script_path)
repository <- dirname(dirname(dirname(verification_dir)))

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Fast-grid verification requires the development package 'pkgload'.")
}
pkgload::load_all(repository, quiet = TRUE)

arguments <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% arguments
unknown <- setdiff(arguments, "--quick")
if (length(unknown) > 0L) {
  stop("Unknown option: ", unknown[[1L]])
}

accuracy_sizes <- if (quick) c(100L, 1000L) else c(100L, 1000L, 5000L)
performance_sizes <- if (quick) c(100L, 1000L) else c(100L, 1000L, 10000L)
grid_n <- 1000L

translation_kernels <- c(
  "gaussian", "epanechnikov", "rectangular", "triangular",
  "biweight", "cosine", "optcosine"
)

translation_scenarios <- list(
  normal = list(
    parametric_start = "normal",
    sample = function(n) stats::rnorm(n)
  ),
  skewed = list(
    parametric_start = "gamma",
    sample = function(n) stats::rgamma(n, shape = 2, rate = 2)
  ),
  bimodal = list(
    parametric_start = "normal",
    sample = function(n) {
      group <- stats::rbinom(n, 1, 0.45)
      stats::rnorm(
        n,
        mean = ifelse(group == 1, -1.8, 1.5),
        sd = ifelse(group == 1, 0.45, 0.7)
      )
    }
  )
)

copula_scenarios <- list(
  interior = function(n) stats::rbeta(n, 2, 5),
  boundary_heavy = function(n) stats::rbeta(n, 0.7, 0.7)
)

trapezoid <- function(x, y) {
  if (length(x) < 2L) {
    return(0)
  }
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

accuracy_metrics <- function(x, fast, exact, rectangular = FALSE) {
  comparable <- is.finite(fast) & is.finite(exact)
  finite_ok <- all(is.finite(fast[is.finite(exact)]))
  if (!any(comparable)) {
    return(data.frame(
      finite = FALSE,
      nonnegative = FALSE,
      scaled_rmse = Inf,
      scaled_max_error = Inf,
      max_central_relative_error = Inf,
      relative_integrated_absolute_error = Inf,
      relative_mass_error = Inf,
      passed = FALSE
    ))
  }

  x <- x[comparable]
  fast <- fast[comparable]
  exact <- exact[comparable]
  peak <- max(exact)
  scale <- max(peak, .Machine$double.eps)
  error <- fast - exact
  central <- exact >= 0.01 * scale
  exact_mass <- trapezoid(x, exact)

  limits <- if (rectangular) {
    c(rmse = 0.01, maximum = 0.04, relative = 0.05, integrated = 0.01)
  } else {
    c(rmse = 0.002, maximum = 0.01, relative = 0.02, integrated = 0.003)
  }

  metrics <- list(
    finite = finite_ok,
    nonnegative = all(fast >= -1e-12 * scale),
    scaled_rmse = sqrt(mean(error^2)) / scale,
    scaled_max_error = max(abs(error)) / scale,
    max_central_relative_error = if (any(central)) {
      max(abs(error[central]) / exact[central])
    } else {
      0
    },
    relative_integrated_absolute_error =
      trapezoid(x, abs(error)) / max(abs(exact_mass), .Machine$double.eps),
    relative_mass_error =
      abs(trapezoid(x, fast) - exact_mass) /
        max(abs(exact_mass), .Machine$double.eps)
  )
  metrics$passed <- with(metrics,
    finite &&
      nonnegative &&
      scaled_rmse <= limits[["rmse"]] &&
      scaled_max_error <= limits[["maximum"]] &&
      max_central_relative_error <= limits[["relative"]] &&
      relative_integrated_absolute_error <= limits[["integrated"]] &&
      relative_mass_error <= limits[["integrated"]]
  )
  as.data.frame(metrics)
}

translation_accuracy_case <- function(kernel, scenario_name, start, n,
                                      bandwidth_ratio) {
  scenario <- translation_scenarios[[scenario_name]]
  seed <- 20260720L + match(kernel, translation_kernels) * 100000L +
    match(scenario_name, names(translation_scenarios)) * 10000L +
    match(start, c("uniform", scenario$parametric_start)) * 1000L +
    n + round(100 * bandwidth_ratio)
  set.seed(seed)
  x <- scenario$sample(n)
  target_h <- bandwidth_ratio * stats::sd(x)
  kernel_sd <- kdensity:::get_kernel(kernel)$sd
  fit <- tryCatch(
    kdensity::kdensity(
      x,
      kernel = kernel,
      start = start,
      bw = target_h / kernel_sd,
      normalized = identical(start, "uniform")
    ),
    error = function(error) {
      stop(
        "Could not construct accuracy case kernel=", kernel,
        ", scenario=", scenario_name,
        ", start=", start,
        ", n=", n,
        ", bandwidth_ratio=", bandwidth_ratio,
        ": ", conditionMessage(error),
        call. = FALSE
      )
    }
  )
  fast <- kdensity::kdensity_grid(fit, n = grid_n, method = "fft")
  exact <- fit(fast$x)
  metrics <- accuracy_metrics(
    fast$x,
    fast$y,
    exact,
    rectangular = kernel == "rectangular"
  )

  cbind(
    data.frame(
      family = "translation",
      kernel = kernel,
      promotion_case = kernel == "gaussian",
      scenario = scenario_name,
      start = start,
      sample_size = n,
      bandwidth = fit$h,
      grid_points = grid_n,
      fft_points = fast$fft_n
    ),
    metrics
  )
}

copula_accuracy_case <- function(scenario_name, start, n, bandwidth) {
  seed <- 20260720L +
    match(scenario_name, names(copula_scenarios)) * 10000L +
    match(start, c("uniform", "beta")) * 1000L +
    n + round(1000 * bandwidth)
  set.seed(seed)
  x <- copula_scenarios[[scenario_name]](n)
  fit <- kdensity::kdensity(
    x,
    kernel = "gcopula",
    start = start,
    bw = bandwidth
  )
  fast <- kdensity::kdensity_grid(fit, n = grid_n, method = "fft")
  exact <- fit(fast$x)
  native <- accuracy_metrics(fast$x, fast$y, exact)

  linear_fast <- kdensity::kdensity_grid(
    fit,
    n = max(2048L, 2L * grid_n),
    method = "fft"
  )
  unique_grid <- !duplicated(linear_fast$x)
  requested <- seq(
    min(linear_fast$x),
    max(linear_fast$x),
    length.out = grid_n
  )
  interpolated <- stats::approx(
    linear_fast$x[unique_grid],
    linear_fast$y[unique_grid],
    xout = requested,
    rule = 2
  )$y
  linear <- accuracy_metrics(requested, interpolated, fit(requested))

  native_metadata <- data.frame(
    family = "copula",
    kernel = "gcopula",
    promotion_case = TRUE,
    scenario = scenario_name,
    start = start,
    sample_size = n,
    bandwidth = fit$h,
    grid_points = grid_n,
    fft_points = fast$fft_n
  )
  linear_metadata <- native_metadata
  linear_metadata$fft_points <- linear_fast$fft_n
  rbind(
    cbind(native_metadata, grid_layout = "native_probit", native),
    cbind(linear_metadata, grid_layout = "requested_linear", linear)
  )
}

translation_cases <- expand.grid(
  kernel = translation_kernels,
  scenario = names(translation_scenarios),
  sample_size = accuracy_sizes,
  bandwidth_ratio = c(0.05, 0.2, 0.6),
  stringsAsFactors = FALSE
)

translation_accuracy <- do.call(rbind, lapply(
  seq_len(nrow(translation_cases)),
  function(index) {
    case <- translation_cases[index, ]
    scenario <- translation_scenarios[[case$scenario]]
    starts <- unique(c("uniform", scenario$parametric_start))
    do.call(rbind, lapply(starts, function(start) {
      translation_accuracy_case(
        case$kernel,
        case$scenario,
        start,
        case$sample_size,
        case$bandwidth_ratio
      )
    }))
  }
))
translation_accuracy$grid_layout <- "linear"

copula_cases <- expand.grid(
  scenario = names(copula_scenarios),
  start = c("uniform", "beta"),
  sample_size = accuracy_sizes,
  bandwidth = c(0.04, 0.12, 0.35),
  stringsAsFactors = FALSE
)
copula_accuracy <- do.call(rbind, lapply(
  seq_len(nrow(copula_cases)),
  function(index) {
    case <- copula_cases[index, ]
    copula_accuracy_case(
      case$scenario,
      case$start,
      case$sample_size,
      case$bandwidth
    )
  }
))

accuracy <- rbind(
  translation_accuracy[, names(copula_accuracy)],
  copula_accuracy
)
row.names(accuracy) <- NULL

adaptive_seconds <- function(fun, minimum_seconds = 0.1, maximum_calls = 512L) {
  calls <- 1L
  repeat {
    elapsed <- unname(system.time(for (index in seq_len(calls)) fun())[["elapsed"]])
    if (elapsed >= minimum_seconds || calls >= maximum_calls) {
      return(c(seconds = elapsed / calls, calls = calls))
    }
    calls <- calls * 2L
  }
}

performance_case <- function(kernel, start, n) {
  set.seed(
    20260720L + match(kernel, c(translation_kernels, "gcopula")) * 100000L +
      match(start, c("uniform", "normal", "beta")) * 10000L + n
  )
  if (kernel == "gcopula") {
    x <- stats::rbeta(n, 2, 5)
    bandwidth <- 0.12
  } else {
    x <- stats::rnorm(n)
    bandwidth <- 0.2 * stats::sd(x) / kdensity:::get_kernel(kernel)$sd
  }
  fit <- kdensity::kdensity(
    x,
    kernel = kernel,
    start = start,
    bw = bandwidth
  )
  fast <- kdensity::kdensity_grid(fit, n = grid_n, method = "fft")
  evaluation_grid <- fast$x

  invisible(fit(evaluation_grid))
  invisible(kdensity::kdensity_grid(fit, n = grid_n, method = "fft"))
  exact_timing <- adaptive_seconds(function() fit(evaluation_grid))
  fft_timing <- adaptive_seconds(function() {
    kdensity::kdensity_grid(fit, n = grid_n, method = "fft")
  })

  exact_seconds <- exact_timing[["seconds"]]
  fft_seconds <- fft_timing[["seconds"]]
  speedup <- exact_seconds / fft_seconds
  promotion_case <- kernel %in% c("gaussian", "gcopula")
  ratio_pass <- if (!promotion_case) {
    TRUE
  } else if (n < 1000L) {
    fft_seconds < 0.01
  } else {
    speedup >= 2
  }

  data.frame(
    kernel = kernel,
    start = start,
    promotion_case = promotion_case,
    sample_size = n,
    grid_points = grid_n,
    exact_seconds = exact_seconds,
    fft_seconds = fft_seconds,
    speedup = speedup,
    exact_calls = exact_timing[["calls"]],
    fft_calls = fft_timing[["calls"]],
    passed = ratio_pass
  )
}

performance_definitions <- rbind(
  expand.grid(
    kernel = translation_kernels,
    start = "uniform",
    sample_size = performance_sizes,
    stringsAsFactors = FALSE
  ),
  expand.grid(
    kernel = "gaussian",
    start = "normal",
    sample_size = performance_sizes,
    stringsAsFactors = FALSE
  ),
  expand.grid(
    kernel = "gcopula",
    start = c("uniform", "beta"),
    sample_size = performance_sizes,
    stringsAsFactors = FALSE
  )
)
performance <- do.call(rbind, lapply(
  seq_len(nrow(performance_definitions)),
  function(index) {
    case <- performance_definitions[index, ]
    performance_case(case$kernel, case$start, case$sample_size)
  }
))
row.names(performance) <- NULL

accuracy_pass <- all(accuracy$passed[accuracy$promotion_case])
performance_pass <- all(performance$passed[performance$promotion_case])
promotion_pass <- accuracy_pass && performance_pass

utils::write.csv(
  accuracy,
  file.path(verification_dir, "accuracy.csv"),
  row.names = FALSE
)
utils::write.csv(
  performance,
  file.path(verification_dir, "performance.csv"),
  row.names = FALSE
)

worst <- function(column, data) {
  index <- which.max(data[[column]])
  paste0(
    format(data[[column]][[index]], digits = 4),
    " (", data$kernel[[index]], ", ", data$scenario[[index]],
    ", start=", data$start[[index]], ", n=", data$sample_size[[index]],
    ", h=", format(data$bandwidth[[index]], digits = 4), ")"
  )
}

promotion_accuracy <- accuracy[accuracy$promotion_case, ]
characterization_accuracy <- accuracy[!accuracy$promotion_case, ]

timing_rows <- vapply(seq_len(nrow(performance)), function(index) {
  row <- performance[index, ]
  paste0(
    "| ", row$kernel, " | ", row$start, " | ", row$sample_size,
    " | ", format(1000 * row$exact_seconds, digits = 4),
    " | ", format(1000 * row$fft_seconds, digits = 4),
    " | ", format(row$speedup, digits = 3),
    " | ", if (row$passed) "yes" else "no", " |"
  )
}, character(1))

report <- c(
  "# Fast-grid verification report",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Decision",
  "",
  paste0(
    if (promotion_pass) "**PASS**" else "**FAIL**",
    ": accuracy ", if (accuracy_pass) "passed" else "failed",
    " and performance ", if (performance_pass) "passed." else "failed."
  ),
  "",
  paste0("- Promotion accuracy cases: ", sum(accuracy$promotion_case)),
  paste0(
    "- Promotion accuracy failures: ",
    sum(!accuracy$passed & accuracy$promotion_case)
  ),
  paste0("- Characterization-only accuracy cases: ", sum(!accuracy$promotion_case)),
  paste0("- Promotion performance cases: ", sum(performance$promotion_case)),
  paste0(
    "- Promotion performance failures: ",
    sum(!performance$passed & performance$promotion_case)
  ),
  "",
  "## Worst promotion accuracy metrics",
  "",
  paste0("- Peak-scaled RMSE: ", worst("scaled_rmse", promotion_accuracy)),
  paste0(
    "- Peak-scaled maximum error: ",
    worst("scaled_max_error", promotion_accuracy)
  ),
  paste0(
    "- Maximum central relative error: ",
    worst("max_central_relative_error", promotion_accuracy)
  ),
  paste0(
    "- Relative integrated absolute error: ",
    worst("relative_integrated_absolute_error", promotion_accuracy)
  ),
  paste0(
    "- Relative mass error: ",
    worst("relative_mass_error", promotion_accuracy)
  ),
  "",
  "## Explicit compact-FFT characterization",
  "",
  paste0(
    "- Worst peak-scaled RMSE: ",
    worst("scaled_rmse", characterization_accuracy)
  ),
  paste0(
    "- Worst relative integrated absolute error: ",
    worst("relative_integrated_absolute_error", characterization_accuracy)
  ),
  "",
  "These cases remain explicit opt-ins; automatic plotting uses their exact",
  "sorted evaluators.",
  "",
  "## Warmed evaluation timings",
  "",
  "| Kernel | Start | n | Exact ms | FFT ms | Speedup | Accepted |",
  "|---|---|---:|---:|---:|---:|:---:|",
  timing_rows,
  "",
  "Full case-level results are in `accuracy.csv` and `performance.csv`."
)
writeLines(report, file.path(verification_dir, "report.md"))

cat(paste(report[1:24], collapse = "\n"), "\n")
if (!promotion_pass) {
  quit(status = 1L)
}
