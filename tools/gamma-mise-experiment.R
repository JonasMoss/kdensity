#!/usr/bin/env Rscript

# Compare gamma-kernel bandwidth selectors by Monte Carlo integrated squared
# error. Run from anywhere with, for example:
#
# Rscript tools/gamma-mise-experiment.R \
#   --reps=20 --n=200 --degree=4 --points=600

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1]]))
script_directory <- dirname(script_path)
repository <- dirname(script_directory)

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("The experiment requires the development package 'pkgload'.")
}
pkgload::load_all(repository, quiet = TRUE)
source(file.path(script_directory, "gamma-start-amise.R"))
source(file.path(script_directory, "mise-oracle-helpers.R"))

parse_options <- function(arguments) {
  options <- list(
    reps = 20L,
    n = 200L,
    degree = 4L,
    points = 600L,
    oracle_points = 0L,
    bootstrap_reps = 199L,
    transition_power = 1,
    seed = 20260719L,
    scenarios = c(
      "gamma",
      "gamma_boundary",
      "gamma_mixture",
      "gamma_close_mixture",
      "lognormal",
      "weibull"
    ),
    methods = c(
      "parametric_gamma",
      "uniform_GR",
      "uniform_GR_fallback",
      "uniform_UCV",
      "gamma_UCV",
      "gamma_Laguerre2",
      "gamma_LaguerreK",
      "gamma_LaguerreK_threshold",
      "gamma_LaguerreK_soft",
      "gamma_LaguerreK_fitted_soft",
      "gamma_LaguerreK_adaptive",
      "gamma_LaguerreK_global",
      "gamma_LaguerreK_conservative",
      "gamma_LaguerreK_bootstrap",
      "gamma_LaguerreK_bootstrap_smooth",
      "gamma_LaguerreK_bic",
      "gamma_LaguerreK_aic",
      "gamma_LaguerreK_garrote",
      "gamma_LaguerreK_raw"
    ),
    output = file.path(repository, "tools", "results")
  )

  for (argument in arguments) {
    pieces <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1]]
    if (length(pieces) != 2L || !pieces[[1]] %in% names(options)) {
      stop("Unknown option: ", argument)
    }
    name <- pieces[[1]]
    value <- pieces[[2]]
    if (name %in% c(
      "reps", "n", "degree", "points", "oracle_points", "bootstrap_reps",
      "seed"
    )) {
      options[[name]] <- as.integer(value)
    } else if (name == "transition_power") {
      options[[name]] <- as.numeric(value)
    } else if (name %in% c("scenarios", "methods")) {
      options[[name]] <- strsplit(value, ",", fixed = TRUE)[[1]]
    } else {
      options[[name]] <- value
    }
  }

  if (options$reps < 1L || options$n < 10L || options$points < 100L) {
    stop("Use reps >= 1, n >= 10, and points >= 100.")
  }
  if (options$degree < 2L) {
    stop("Use degree >= 2.")
  }
  if (options$oracle_points == 1L || options$oracle_points < 0L) {
    stop("Use oracle_points = 0 to disable or oracle_points >= 2.")
  }
  if (options$bootstrap_reps < 19L) {
    stop("Use bootstrap_reps >= 19.")
  }
  if (!is.finite(options$transition_power) || options$transition_power <= 0) {
    stop("Use transition_power > 0.")
  }
  if (any(grepl("LaguerreK_bootstrap", options$methods, fixed = TRUE))) {
    minimum_bootstrap_reps <- ceiling(
      1 / conservative_series_alpha(options$n)
    ) - 1L
    if (options$bootstrap_reps < minimum_bootstrap_reps) {
      stop(
        "The shrinking test level is unattainable with bootstrap_reps = ",
        options$bootstrap_reps,
        "; use at least ",
        minimum_bootstrap_reps,
        " for n = ",
        options$n,
        "."
      )
    }
  }
  options
}

mixture_sample <- function(n, probability, first, second) {
  use_first <- stats::runif(n) < probability
  out <- numeric(n)
  out[use_first] <- first(sum(use_first))
  out[!use_first] <- second(sum(!use_first))
  out
}

scenario_catalog <- list(
  gamma = list(
    sample = function(n) stats::rgamma(n, shape = 3, rate = 1.5),
    density = function(x) stats::dgamma(x, shape = 3, rate = 1.5),
    upper = stats::qgamma(0.9999, shape = 3, rate = 1.5)
  ),
  gamma_boundary = list(
    sample = function(n) stats::rgamma(n, shape = 0.8, rate = 1),
    density = function(x) stats::dgamma(x, shape = 0.8, rate = 1),
    upper = stats::qgamma(0.9999, shape = 0.8, rate = 1)
  ),
  gamma_mixture = list(
    sample = function(n) {
      mixture_sample(
        n,
        probability = 0.65,
        first = function(size) stats::rgamma(size, shape = 2, rate = 1.2),
        second = function(size) stats::rgamma(size, shape = 12, rate = 2)
      )
    },
    density = function(x) {
      0.65 * stats::dgamma(x, shape = 2, rate = 1.2) +
        0.35 * stats::dgamma(x, shape = 12, rate = 2)
    },
    upper = max(
      stats::qgamma(0.9999, shape = 2, rate = 1.2),
      stats::qgamma(0.9999, shape = 12, rate = 2)
    )
  ),
  gamma_close_mixture = list(
    sample = function(n) {
      mixture_sample(
        n,
        probability = 0.7,
        first = function(size) stats::rgamma(size, shape = 3, rate = 1.5),
        second = function(size) stats::rgamma(size, shape = 6, rate = 3)
      )
    },
    density = function(x) {
      0.7 * stats::dgamma(x, shape = 3, rate = 1.5) +
        0.3 * stats::dgamma(x, shape = 6, rate = 3)
    },
    upper = max(
      stats::qgamma(0.9999, shape = 3, rate = 1.5),
      stats::qgamma(0.9999, shape = 6, rate = 3)
    )
  ),
  lognormal = list(
    sample = function(n) stats::rlnorm(n, meanlog = 0.4, sdlog = 0.65),
    density = function(x) stats::dlnorm(x, meanlog = 0.4, sdlog = 0.65),
    upper = stats::qlnorm(0.9999, meanlog = 0.4, sdlog = 0.65)
  ),
  weibull = list(
    sample = function(n) stats::rweibull(n, shape = 1.4, scale = 2),
    density = function(x) stats::dweibull(x, shape = 1.4, scale = 2),
    upper = stats::qweibull(0.9999, shape = 1.4, scale = 2)
  )
)

trapezoid_integral <- function(x, y) {
  quadrature_weights <- attr(x, "quadrature_weights")
  if (!is.null(quadrature_weights)) {
    return(sum(quadrature_weights * y))
  }
  sum(diff(x) * (y[-length(y)] + y[-1L]) / 2)
}

make_evaluation_grid <- function(scenario, points) {
  u <- (seq_len(points) - 0.5) / points
  grid <- scenario$upper * u^2
  attr(grid, "quadrature_weights") <- 2 * scenario$upper * u / points
  grid
}

capture_warnings <- function(expression) {
  warnings <- character()
  value <- withCallingHandlers(
    expression,
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warnings)
}

fit_method <- function(method, x, degree, bootstrap_reps, transition_power) {
  selector <- NULL
  started <- proc.time()[["elapsed"]]

  result <- tryCatch(
    capture_warnings({
      fit <- switch(method,
        parametric_gamma = kdensity::kdensity(
          x, kernel = "gamma", start = "gamma", bw = Inf
        ),
        uniform_GR = {
          bandwidth <- kdensity:::compute_gr_bandwidth(x)
          kdensity::kdensity(
            x, kernel = "gamma", start = "uniform", bw = bandwidth
          )
        },
        uniform_GR_fallback = kdensity::kdensity(
          x, kernel = "gamma", start = "uniform", bw = "GR"
        ),
        uniform_UCV = kdensity::kdensity(
          x, kernel = "gamma", start = "uniform", bw = "ucv"
        ),
        gamma_UCV = kdensity::kdensity(
          x, kernel = "gamma", start = "gamma", bw = "ucv"
        ),
        gamma_Laguerre2 = {
          selector <- laguerre_gamma_bandwidth(
            x, degree = 2L, normalized = TRUE
          )
          kdensity::kdensity(
            x, kernel = "gamma", start = "gamma",
            bw = selector$bandwidth
          )
        },
        gamma_LaguerreK = {
          selector <- laguerre_gamma_bandwidth(
            x, degree = degree, normalized = TRUE
          )
          kdensity::kdensity(
            x, kernel = "gamma", start = "gamma",
            bw = selector$bandwidth
          )
        },
        gamma_LaguerreK_threshold = {
          selector <- laguerre_gamma_bandwidth(
            x, degree = degree, normalized = TRUE, threshold = 1
          )
          kdensity::kdensity(
            x, kernel = "gamma", start = "gamma",
            bw = selector$bandwidth
          )
        },
        gamma_LaguerreK_soft = {
          selector <- laguerre_gamma_bandwidth(
            x,
            degree = degree,
            normalized = TRUE,
            threshold = 1,
            shrinkage = "soft"
          )
          kdensity::kdensity(
            x, kernel = "gamma", start = "gamma",
            bw = selector$bandwidth
          )
        },
        gamma_LaguerreK_fitted_soft = {
          selector <- laguerre_gamma_bandwidth(
            x,
            degree = degree,
            normalized = TRUE,
            threshold = 1,
            shrinkage = "soft",
            threshold_scale = "fitted"
          )
          kdensity::kdensity(
            x, kernel = "gamma", start = "gamma",
            bw = selector$bandwidth
          )
        },
        gamma_LaguerreK_adaptive = {
          selector <- laguerre_gamma_bandwidth(
            x,
            degree = degree,
            normalized = TRUE,
            threshold = sqrt(2 * log(max(degree - 1L, 2L))),
            shrinkage = "hard",
            threshold_scale = "fitted"
          )
          kdensity::kdensity(
            x, kernel = "gamma", start = "gamma",
            bw = selector$bandwidth
          )
        },
        gamma_LaguerreK_global = {
          selector <- laguerre_gamma_bandwidth(
            x,
            degree = degree,
            normalized = TRUE,
            threshold_scale = "fitted"
          )
          selector <- apply_global_series_test(selector, length(x))
          kdensity::kdensity(
            x, kernel = "gamma", start = "gamma",
            bw = selector$bandwidth
          )
        },
        gamma_LaguerreK_conservative = {
          selector <- laguerre_gamma_bandwidth(
            x,
            degree = degree,
            normalized = TRUE,
            threshold_scale = "fitted"
          )
          selector <- apply_global_series_test(
            selector,
            length(x),
            alpha = conservative_series_alpha(length(x))
          )
          kdensity::kdensity(
            x, kernel = "gamma", start = "gamma",
            bw = selector$bandwidth
          )
        },
        gamma_LaguerreK_bootstrap = {
          selector <- laguerre_gamma_bandwidth(
            x,
            degree = degree,
            normalized = TRUE,
            threshold_scale = "fitted"
          )
          shape <- selector$shape
          rate <- selector$rate
          selector <- apply_bootstrap_global_series_test(
            selector,
            length(x),
            simulate = function(size) {
              stats::rgamma(size, shape = shape, rate = rate)
            },
            refit = function(bootstrap_x) {
              laguerre_gamma_bandwidth(
                bootstrap_x,
                degree = degree,
                normalized = TRUE,
                threshold_scale = "fitted"
              )
            },
            replicates = bootstrap_reps
          )
          kdensity::kdensity(
            x, kernel = "gamma", start = "gamma",
            bw = selector$bandwidth
          )
        },
        gamma_LaguerreK_bootstrap_smooth = {
          selector <- laguerre_gamma_bandwidth(
            x,
            degree = degree,
            normalized = TRUE,
            threshold_scale = "fitted"
          )
          shape <- selector$shape
          rate <- selector$rate
          selector <- apply_bootstrap_global_series_test(
            selector,
            length(x),
            simulate = function(size) {
              stats::rgamma(size, shape = shape, rate = rate)
            },
            refit = function(bootstrap_x) {
              laguerre_gamma_bandwidth(
                bootstrap_x,
                degree = degree,
                normalized = TRUE,
                threshold_scale = "fitted"
              )
            },
            replicates = bootstrap_reps
          )
          selector <- smooth_bootstrap_bandwidth(
            selector,
            power = transition_power
          )
          kdensity::kdensity(
            x, kernel = "gamma", start = "gamma",
            bw = selector$bandwidth
          )
        },
        gamma_LaguerreK_bic = {
          preliminary <- laguerre_gamma_bandwidth(
            x,
            degree = degree,
            normalized = TRUE,
            threshold_scale = "fitted"
          )
          selection <- select_series_degree(preliminary, length(x))
          if (selection$degree == 0L) {
            selector <- preliminary
            selector$bandwidth <- Inf
          } else {
            selector <- laguerre_gamma_bandwidth(
              x,
              degree = selection$degree,
              normalized = TRUE,
              threshold_scale = "fitted"
            )
          }
          selector$selected_degree <- selection$degree
          selector$degree_scores <- selection$scores
          kdensity::kdensity(
            x, kernel = "gamma", start = "gamma",
            bw = selector$bandwidth
          )
        },
        gamma_LaguerreK_aic = {
          preliminary <- laguerre_gamma_bandwidth(
            x,
            degree = degree,
            normalized = TRUE,
            threshold_scale = "fitted"
          )
          selection <- select_series_degree(
            preliminary,
            length(x),
            penalty = 2
          )
          if (selection$degree == 0L) {
            selector <- preliminary
            selector$bandwidth <- Inf
          } else {
            selector <- laguerre_gamma_bandwidth(
              x,
              degree = selection$degree,
              normalized = TRUE,
              threshold_scale = "fitted"
            )
          }
          selector$selected_degree <- selection$degree
          selector$degree_scores <- selection$scores
          kdensity::kdensity(
            x, kernel = "gamma", start = "gamma",
            bw = selector$bandwidth
          )
        },
        gamma_LaguerreK_garrote = {
          selector <- laguerre_gamma_bandwidth(
            x,
            degree = degree,
            normalized = TRUE,
            threshold = 1,
            shrinkage = "garrote"
          )
          kdensity::kdensity(
            x, kernel = "gamma", start = "gamma",
            bw = selector$bandwidth
          )
        },
        gamma_LaguerreK_raw = {
          selector <- laguerre_gamma_bandwidth(
            x, degree = degree, normalized = FALSE
          )
          kdensity::kdensity(
            x, kernel = "gamma", start = "gamma",
            bw = selector$bandwidth
          )
        },
        stop("Unknown method: ", method)
      )
      fit
    }),
    error = function(condition) {
      list(error = conditionMessage(condition))
    }
  )

  elapsed <- proc.time()[["elapsed"]] - started
  if (!is.null(result$error)) {
    return(list(
      fit = NULL,
      elapsed = elapsed,
      warning = "",
      error = result$error,
      selector = selector
    ))
  }

  list(
    fit = result$value,
    elapsed = elapsed,
    warning = paste(result$warnings, collapse = " | "),
    error = "",
    selector = selector
  )
}

one_result <- function(scenario_name, replication, x, method, degree, points,
                       bootstrap_reps, transition_power) {
  scenario <- scenario_catalog[[scenario_name]]
  n <- length(x)
  grid <- make_evaluation_grid(scenario, points)
  truth <- scenario$density(grid)
  fitted <- fit_method(
    method,
    x,
    degree,
    bootstrap_reps,
    transition_power
  )

  ise <- NA_real_
  bandwidth <- NA_real_
  theta2 <- NA_real_
  theta3 <- NA_real_
  theta4 <- NA_real_
  roughness <- NA_real_
  normalization_bias <- NA_real_
  global_p_value <- NA_real_
  global_asymptotic_p_value <- NA_real_
  global_selected <- NA
  bootstrap_successful <- NA_integer_
  bootstrap_failures <- NA_integer_
  bootstrap_exceedances <- NA_integer_
  bootstrap_critical_value <- NA_real_
  smooth_weight <- NA_real_
  selected_degree <- NA_integer_

  if (!is.null(fitted$fit)) {
    estimate <- tryCatch(
      fitted$fit(grid),
      error = function(condition) rep(NA_real_, length(grid))
    )
    if (all(is.finite(estimate))) {
      ise <- trapezoid_integral(grid, (estimate - truth)^2)
    } else if (!nzchar(fitted$error)) {
      fitted$error <- "Density evaluation returned non-finite values."
    }
    bandwidth <- fitted$fit$bw
  }

  if (!is.null(fitted$selector)) {
    selector <- fitted$selector
    theta_value <- function(name) {
      if (name %in% names(selector$theta)) {
        unname(selector$theta[[name]])
      } else {
        NA_real_
      }
    }
    theta2 <- theta_value("theta2")
    theta3 <- theta_value("theta3")
    theta4 <- theta_value("theta4")
    roughness <- selector$roughness
    normalization_bias <- selector$normalization_bias
    if (!is.null(selector$global_p_value)) {
      global_p_value <- selector$global_p_value
      global_selected <- selector$global_selected
    }
    if (!is.null(selector$global_asymptotic_p_value)) {
      global_asymptotic_p_value <- selector$global_asymptotic_p_value
    }
    if (!is.null(selector$bootstrap_replicates)) {
      bootstrap_successful <- selector$bootstrap_successful
      bootstrap_failures <- selector$bootstrap_failures
      bootstrap_exceedances <- selector$bootstrap_exceedances
      bootstrap_critical_value <- selector$bootstrap_critical_value
    }
    if (!is.null(selector$smooth_weight)) {
      smooth_weight <- selector$smooth_weight
    }
    if (!is.null(selector$selected_degree)) {
      selected_degree <- selector$selected_degree
    }
  }

  data.frame(
    scenario = scenario_name,
    replication = replication,
    n = n,
    method = method,
    degree = if (method == "gamma_Laguerre2") {
      2L
    } else if (grepl("LaguerreK", method, fixed = TRUE)) {
      degree
    } else {
      NA_integer_
    },
    bandwidth = bandwidth,
    ise = ise,
    fit_seconds = fitted$elapsed,
    theta2 = theta2,
    theta3 = theta3,
    theta4 = theta4,
    roughness = roughness,
    normalization_bias = normalization_bias,
    global_p_value = global_p_value,
    global_asymptotic_p_value = global_asymptotic_p_value,
    global_selected = global_selected,
    bootstrap_successful = bootstrap_successful,
    bootstrap_failures = bootstrap_failures,
    bootstrap_exceedances = bootstrap_exceedances,
    bootstrap_critical_value = bootstrap_critical_value,
    smooth_weight = smooth_weight,
    selected_degree = selected_degree,
    warning = fitted$warning,
    error = fitted$error,
    stringsAsFactors = FALSE
  )
}

summarize_results <- function(results) {
  groups <- split(
    results,
    interaction(results$scenario, results$method, drop = TRUE)
  )
  summaries <- lapply(groups, function(group) {
    successful <- is.finite(group$ise)
    ise <- group$ise[successful]
    data.frame(
      scenario = group$scenario[[1]],
      method = group$method[[1]],
      n = group$n[[1]],
      repetitions = nrow(group),
      successful = sum(successful),
      failures = sum(!successful),
      mise = if (length(ise)) mean(ise) else NA_real_,
      mise_se = if (length(ise) > 1L) {
        stats::sd(ise) / sqrt(length(ise))
      } else {
        NA_real_
      },
      median_bandwidth = if (any(is.finite(group$bandwidth))) {
        stats::median(group$bandwidth[is.finite(group$bandwidth)])
      } else {
        NA_real_
      },
      mean_fit_seconds = mean(group$fit_seconds),
      global_selections = sum(group$global_selected, na.rm = TRUE),
      bootstrap_failures = sum(group$bootstrap_failures, na.rm = TRUE),
      mean_smooth_weight = if (any(is.finite(group$smooth_weight))) {
        mean(group$smooth_weight, na.rm = TRUE)
      } else {
        NA_real_
      },
      median_selected_degree = if (any(!is.na(group$selected_degree))) {
        stats::median(group$selected_degree, na.rm = TRUE)
      } else {
        NA_real_
      },
      fallbacks = sum(grepl("Falling back", group$warning, fixed = TRUE)),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, summaries)
  summary[order(summary$scenario, summary$mise), ]
}

options <- parse_options(commandArgs(trailingOnly = TRUE))
unknown_scenarios <- setdiff(options$scenarios, names(scenario_catalog))
if (length(unknown_scenarios)) {
  stop("Unknown scenarios: ", paste(unknown_scenarios, collapse = ", "))
}

set.seed(options$seed)
dir.create(options$output, recursive = TRUE, showWarnings = FALSE)
total <- options$reps * length(options$scenarios) * length(options$methods)
results <- vector("list", total)
index <- 0L
oracle_results <- list()
oracle_index <- 0L

for (scenario_name in options$scenarios) {
  for (replication in seq_len(options$reps)) {
    # Draw once per scenario/replication so every method sees the same sample.
    scenario <- scenario_catalog[[scenario_name]]
    sample_x <- scenario$sample(options$n)
    message(
      sprintf(
        "[%s] replication %d/%d",
        scenario_name, replication, options$reps
      )
    )

    for (method in options$methods) {
      index <- index + 1L
      results[[index]] <- one_result(
        scenario_name,
        replication,
        sample_x,
        method,
        options$degree,
        options$points,
        options$bootstrap_reps,
        options$transition_power
      )
    }

    if (options$oracle_points >= 2L) {
      grid <- make_evaluation_grid(scenario, options$points)
      truth <- scenario$density(grid)
      finite_bandwidths <- log_bandwidth_grid(
        lower = scenario$upper * 1e-5,
        upper = scenario$upper * 0.2,
        points = options$oracle_points
      )

      oracle_index <- oracle_index + 1L
      oracle_results[[oracle_index]] <- evaluate_oracle_grid(
        sample_x,
        evaluation_grid = grid,
        truth = truth,
        bandwidths = finite_bandwidths,
        fit_function = function(bandwidth) {
          kdensity::kdensity(
            sample_x,
            kernel = "gamma",
            start = "uniform",
            bw = bandwidth
          )
        },
        scenario = scenario_name,
        replication = replication,
        estimator = "uniform_start"
      )

      oracle_index <- oracle_index + 1L
      oracle_results[[oracle_index]] <- evaluate_oracle_grid(
        sample_x,
        evaluation_grid = grid,
        truth = truth,
        bandwidths = c(finite_bandwidths, Inf),
        fit_function = function(bandwidth) {
          kdensity::kdensity(
            sample_x,
            kernel = "gamma",
            start = "gamma",
            bw = bandwidth
          )
        },
        scenario = scenario_name,
        replication = replication,
        estimator = "gamma_start"
      )
    }
  }
}

results <- do.call(rbind, results)
summary <- summarize_results(results)
stem <- sprintf(
  "gamma-mise-n%d-r%d-d%d-seed%d",
  options$n, options$reps, options$degree, options$seed
)
results_path <- file.path(options$output, paste0(stem, ".csv"))
summary_path <- file.path(options$output, paste0(stem, "-summary.csv"))
utils::write.csv(results, results_path, row.names = FALSE)
utils::write.csv(summary, summary_path, row.names = FALSE)

print(summary, row.names = FALSE)
message("Wrote: ", results_path)
message("Wrote: ", summary_path)

if (length(oracle_results)) {
  oracle_results <- do.call(rbind, oracle_results)
  oracle_summary <- summarize_oracle_grid(oracle_results)
  oracle_optima <- oracle_grid_optima(oracle_summary)
  oracle_path <- file.path(options$output, paste0(stem, "-oracle.csv"))
  oracle_summary_path <- file.path(
    options$output,
    paste0(stem, "-oracle-summary.csv")
  )
  oracle_optima_path <- file.path(
    options$output,
    paste0(stem, "-oracle-optima.csv")
  )
  utils::write.csv(oracle_results, oracle_path, row.names = FALSE)
  utils::write.csv(oracle_summary, oracle_summary_path, row.names = FALSE)
  utils::write.csv(oracle_optima, oracle_optima_path, row.names = FALSE)

  message("Oracle-MISE grid optima:")
  print(oracle_optima, row.names = FALSE)
  message("Wrote: ", oracle_path)
  message("Wrote: ", oracle_summary_path)
  message("Wrote: ", oracle_optima_path)
}
