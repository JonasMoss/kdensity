#!/usr/bin/env Rscript

# Compare beta-kernel bandwidth selectors by Monte Carlo integrated squared
# error. Run from anywhere with, for example:
#
# Rscript tools/beta-mise-experiment.R \
#   --reps=20 --n=200 --degree=4 --points=600

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1]]))
script_directory <- dirname(script_path)
repository <- dirname(script_directory)

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("The experiment requires the development package 'pkgload'.")
}
pkgload::load_all(repository, quiet = TRUE)
source(file.path(script_directory, "beta-start-amise.R"))
source(file.path(script_directory, "mise-oracle-helpers.R"))

parse_options <- function(arguments) {
  options <- list(
    reps = 20L,
    n = 200L,
    degree = 4L,
    points = 600L,
    oracle_points = 0L,
    bootstrap_reps = 199L,
    transition_power = 4,
    seed = 20260719L,
    scenarios = c(
      "beta",
      "beta_boundary",
      "beta_mixture",
      "beta_interior_mixture",
      "kumaraswamy",
      "logistic_normal"
    ),
    methods = c(
      "parametric_beta",
      "uniform_HS",
      "uniform_HS_normalized",
      "uniform_HS_fallback",
      "uniform_UCV",
      "beta_UCV",
      "beta_Jacobi2",
      "beta_JacobiK",
      "beta_JacobiK_threshold",
      "beta_JacobiK_soft",
      "beta_JacobiK_fitted_soft",
      "beta_JacobiK_adaptive",
      "beta_JacobiK_global",
      "beta_JacobiK_guarded",
      "beta_JacobiK_bootstrap_guarded",
      "beta_JacobiK_bootstrap_smooth_guarded",
      "beta_JacobiK_bic",
      "beta_JacobiK_aic",
      "beta_JacobiK_garrote",
      "beta_JacobiK_raw"
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
  if (any(grepl("JacobiK_bootstrap", options$methods, fixed = TRUE))) {
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
  beta = list(
    sample = function(n) stats::rbeta(n, shape1 = 3, shape2 = 5),
    density = function(x) stats::dbeta(x, shape1 = 3, shape2 = 5)
  ),
  beta_boundary = list(
    sample = function(n) stats::rbeta(n, shape1 = 0.8, shape2 = 0.8),
    density = function(x) stats::dbeta(x, shape1 = 0.8, shape2 = 0.8)
  ),
  beta_mixture = list(
    sample = function(n) {
      mixture_sample(
        n,
        probability = 0.55,
        first = function(size) stats::rbeta(size, shape1 = 2, shape2 = 8),
        second = function(size) stats::rbeta(size, shape1 = 8, shape2 = 2)
      )
    },
    density = function(x) {
      0.55 * stats::dbeta(x, shape1 = 2, shape2 = 8) +
      0.45 * stats::dbeta(x, shape1 = 8, shape2 = 2)
    }
  ),
  beta_interior_mixture = list(
    sample = function(n) {
      mixture_sample(
        n,
        probability = 0.5,
        first = function(size) stats::rbeta(size, shape1 = 8, shape2 = 20),
        second = function(size) stats::rbeta(size, shape1 = 20, shape2 = 8)
      )
    },
    density = function(x) {
      0.5 * stats::dbeta(x, shape1 = 8, shape2 = 20) +
        0.5 * stats::dbeta(x, shape1 = 20, shape2 = 8)
    }
  ),
  kumaraswamy = list(
    sample = function(n) {
      uniform <- stats::runif(n)
      (1 - (1 - uniform)^(1 / 3))^(1 / 1.6)
    },
    density = function(x) {
      1.6 * 3 * x^0.6 * (1 - x^1.6)^2
    }
  ),
  logistic_normal = list(
    sample = function(n) stats::plogis(stats::rnorm(n, mean = 0.4, sd = 1.2)),
    density = function(x) {
      stats::dnorm(stats::qlogis(x), mean = 0.4, sd = 1.2) /
        (x * (1 - x))
    }
  )
)

trapezoid_integral <- function(x, y) {
  quadrature_weights <- attr(x, "quadrature_weights")
  if (!is.null(quadrature_weights)) {
    return(sum(quadrature_weights * y))
  }
  sum(diff(x) * (y[-length(y)] + y[-1L]) / 2)
}

make_evaluation_grid <- function(points) {
  # Midpoint quadrature after x = sin(pi u / 2)^2 resolves both boundaries
  # without evaluating at zero or one. Retaining the Jacobian avoids the
  # unstable first trapezoid that otherwise occurs for singular densities.
  u <- (seq_len(points) - 0.5) / points
  grid <- sin(pi * u / 2)^2
  attr(grid, "quadrature_weights") <- (pi / 2) * sin(pi * u) / points
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

estimate_endpoint_shapes <- function(x) {
  # If F(x) is proportional to x^a near zero, the lower order statistics
  # give the usual Pareto-tail estimate of a after reversing the ratio.
  k <- max(5L, floor(sqrt(length(x))))
  lower <- sort(x)
  upper <- sort(1 - x)
  c(
    lower = k / sum(log(lower[k + 1L] / lower[seq_len(k)])),
    upper = k / sum(log(upper[k + 1L] / upper[seq_len(k)]))
  )
}

fit_jacobi_method <- function(x, selector) {
  # Chen's corrected boundary construction has disjoint lower/interior/upper
  # regions only for b < 1/4. The AMISE rule is asymptotic and can exceed that
  # limit when its estimated roughness is very small.
  bandwidth <- selector$bandwidth
  if (!is.infinite(bandwidth)) {
    bandwidth <- min(bandwidth, 0.25 - sqrt(.Machine$double.eps))
  }
  kdensity::kdensity(
    x,
    kernel = "beta",
    start = "beta",
    bw = bandwidth
  )
}

fit_guarded_jacobi_method <- function(x, degree, bootstrap_reps = NULL,
                                      smooth = FALSE,
                                      transition_power = 4) {
  start_fit <- beta_start_fit(x)
  # A fitted shape at or below one means that the beta start is flat or
  # singular at an endpoint. Leave a small finite-sample margin around that
  # boundary, but do not reject benign fitted shapes merely because the older
  # closed-form roughness rule needs shapes > 3/2.
  if (min(start_fit) <= 1.1) {
    endpoint_shapes <- estimate_endpoint_shapes(x)
    boundary_bandwidth <- if (min(endpoint_shapes) <= 1) {
      "UCV"
    } else {
      "HS"
    }
    selector <- list(
      bandwidth = NA_real_,
      degree = degree,
      theta = stats::setNames(
        rep(NA_real_, degree + 1L),
        paste0("theta", 0:degree)
      ),
      roughness = NA_real_,
      normalization_bias = NA_real_,
      shape1 = start_fit[["shape1"]],
      shape2 = start_fit[["shape2"]],
      boundary_guard = TRUE,
      boundary_bandwidth = boundary_bandwidth,
      endpoint_shape_lower = endpoint_shapes[["lower"]],
      endpoint_shape_upper = endpoint_shapes[["upper"]]
    )
    fit <- if (boundary_bandwidth == "UCV") {
      kdensity::kdensity(
        x, kernel = "beta", start = "uniform", bw = "ucv"
      )
    } else {
      bandwidth <- kdensity:::compute_hs_bandwidth(x)
      kdensity::kdensity(
        x, kernel = "beta", start = "uniform", bw = bandwidth
      )
    }
  } else {
    selector <- jacobi_beta_bandwidth(
      x,
      degree = degree,
      normalized = TRUE,
      threshold_scale = "fitted"
    )
    if (is.null(bootstrap_reps)) {
      selector <- apply_global_series_test(
        selector,
        length(x),
        alpha = conservative_series_alpha(length(x))
      )
    } else {
      shape1 <- selector$shape1
      shape2 <- selector$shape2
      selector <- apply_bootstrap_global_series_test(
        selector,
        length(x),
        simulate = function(size) {
          stats::rbeta(size, shape1 = shape1, shape2 = shape2)
        },
        refit = function(bootstrap_x) {
          jacobi_beta_bandwidth(
            bootstrap_x,
            degree = degree,
            normalized = TRUE,
            threshold_scale = "fitted"
          )
        },
        replicates = bootstrap_reps
      )
    }
    selector$boundary_guard <- FALSE
    if (smooth) {
      weight <- bootstrap_transition_weight(
        selector,
        power = transition_power
      )
      selector$smooth_weight <- weight
      selector$transition_power <- transition_power
      selector$bandwidth <- selector$unregularized_bandwidth
      parametric_fit <- kdensity::kdensity(
        x, kernel = "beta", start = "beta", bw = Inf
      )
      nonparametric_fit <- fit_jacobi_method(x, selector)
      fit <- weighted_density_fit(
        parametric_fit,
        nonparametric_fit,
        weight
      )
    } else {
      fit <- fit_jacobi_method(x, selector)
    }
  }
  list(fit = fit, selector = selector)
}

fit_method <- function(method, x, degree, bootstrap_reps, transition_power) {
  selector <- NULL
  started <- proc.time()[["elapsed"]]

  result <- tryCatch(
    capture_warnings({
      fit <- switch(method,
        parametric_beta = kdensity::kdensity(
          x, kernel = "beta", start = "beta", bw = Inf
        ),
        uniform_HS = {
          bandwidth <- kdensity:::compute_hs_bandwidth(x)
          kdensity::kdensity(
            x, kernel = "beta", start = "uniform", bw = bandwidth
          )
        },
        uniform_HS_normalized = {
          selector <- normalized_beta_reference_bandwidth(x)
          kdensity::kdensity(
            x, kernel = "beta", start = "uniform",
            bw = selector$bandwidth
          )
        },
        uniform_HS_fallback = kdensity::kdensity(
          x, kernel = "beta", start = "uniform", bw = "HS"
        ),
        uniform_UCV = kdensity::kdensity(
          x, kernel = "beta", start = "uniform", bw = "ucv"
        ),
        beta_UCV = kdensity::kdensity(
          x, kernel = "beta", start = "beta", bw = "ucv"
        ),
        beta_Jacobi2 = {
          selector <- jacobi_beta_bandwidth(
            x, degree = 2L, normalized = TRUE
          )
          fit_jacobi_method(x, selector)
        },
        beta_JacobiK = {
          selector <- jacobi_beta_bandwidth(
            x, degree = degree, normalized = TRUE
          )
          fit_jacobi_method(x, selector)
        },
        beta_JacobiK_threshold = {
          selector <- jacobi_beta_bandwidth(
            x, degree = degree, normalized = TRUE, threshold = 1
          )
          fit_jacobi_method(x, selector)
        },
        beta_JacobiK_soft = {
          selector <- jacobi_beta_bandwidth(
            x,
            degree = degree,
            normalized = TRUE,
            threshold = 1,
            shrinkage = "soft"
          )
          fit_jacobi_method(x, selector)
        },
        beta_JacobiK_fitted_soft = {
          selector <- jacobi_beta_bandwidth(
            x,
            degree = degree,
            normalized = TRUE,
            threshold = 1,
            shrinkage = "soft",
            threshold_scale = "fitted"
          )
          fit_jacobi_method(x, selector)
        },
        beta_JacobiK_adaptive = {
          selector <- jacobi_beta_bandwidth(
            x,
            degree = degree,
            normalized = TRUE,
            threshold = sqrt(2 * log(max(degree - 1L, 2L))),
            shrinkage = "hard",
            threshold_scale = "fitted"
          )
          fit_jacobi_method(x, selector)
        },
        beta_JacobiK_global = {
          selector <- jacobi_beta_bandwidth(
            x,
            degree = degree,
            normalized = TRUE,
            threshold_scale = "fitted"
          )
          selector <- apply_global_series_test(selector, length(x))
          fit_jacobi_method(x, selector)
        },
        beta_JacobiK_guarded = {
          guarded <- fit_guarded_jacobi_method(x, degree)
          selector <- guarded$selector
          guarded$fit
        },
        beta_JacobiK_bootstrap_guarded = {
          guarded <- fit_guarded_jacobi_method(
            x,
            degree,
            bootstrap_reps = bootstrap_reps
          )
          selector <- guarded$selector
          guarded$fit
        },
        beta_JacobiK_bootstrap_smooth_guarded = {
          guarded <- fit_guarded_jacobi_method(
            x,
            degree,
            bootstrap_reps = bootstrap_reps,
            smooth = TRUE,
            transition_power = transition_power
          )
          selector <- guarded$selector
          guarded$fit
        },
        beta_JacobiK_bic = {
          preliminary <- jacobi_beta_bandwidth(
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
            selector <- jacobi_beta_bandwidth(
              x,
              degree = selection$degree,
              normalized = TRUE,
              threshold_scale = "fitted"
            )
          }
          selector$selected_degree <- selection$degree
          selector$degree_scores <- selection$scores
          fit_jacobi_method(x, selector)
        },
        beta_JacobiK_aic = {
          preliminary <- jacobi_beta_bandwidth(
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
            selector <- jacobi_beta_bandwidth(
              x,
              degree = selection$degree,
              normalized = TRUE,
              threshold_scale = "fitted"
            )
          }
          selector$selected_degree <- selection$degree
          selector$degree_scores <- selection$scores
          fit_jacobi_method(x, selector)
        },
        beta_JacobiK_garrote = {
          selector <- jacobi_beta_bandwidth(
            x,
            degree = degree,
            normalized = TRUE,
            threshold = 1,
            shrinkage = "garrote"
          )
          fit_jacobi_method(x, selector)
        },
        beta_JacobiK_raw = {
          selector <- jacobi_beta_bandwidth(
            x, degree = degree, normalized = FALSE
          )
          fit_jacobi_method(x, selector)
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
  grid <- make_evaluation_grid(points)
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
  raw_bandwidth <- NA_real_
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
  shape1 <- NA_real_
  shape2 <- NA_real_
  capped <- FALSE
  boundary_guard <- FALSE
  boundary_bandwidth <- NA_character_
  endpoint_shape_lower <- NA_real_
  endpoint_shape_upper <- NA_real_

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
    shape1 <- selector$shape1
    shape2 <- selector$shape2
    raw_bandwidth <- selector$bandwidth
    capped <- is.finite(raw_bandwidth) &&
      raw_bandwidth > 0.25 - sqrt(.Machine$double.eps)
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
    if (!is.null(selector$boundary_guard)) {
      boundary_guard <- selector$boundary_guard
    }
    if (!is.null(selector$boundary_bandwidth)) {
      boundary_bandwidth <- selector$boundary_bandwidth
      endpoint_shape_lower <- selector$endpoint_shape_lower
      endpoint_shape_upper <- selector$endpoint_shape_upper
    }
  }

  data.frame(
    scenario = scenario_name,
    replication = replication,
    n = n,
    method = method,
    degree = if (method == "beta_Jacobi2") {
      2L
    } else if (grepl("JacobiK", method, fixed = TRUE)) {
      degree
    } else {
      NA_integer_
    },
    bandwidth = bandwidth,
    raw_bandwidth = raw_bandwidth,
    capped = capped,
    boundary_guard = boundary_guard,
    boundary_bandwidth = boundary_bandwidth,
    endpoint_shape_lower = endpoint_shape_lower,
    endpoint_shape_upper = endpoint_shape_upper,
    ise = ise,
    fit_seconds = fitted$elapsed,
    shape1 = shape1,
    shape2 = shape2,
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
      caps = sum(group$capped),
      boundary_guards = sum(group$boundary_guard),
      boundary_ucv = sum(group$boundary_bandwidth == "UCV", na.rm = TRUE),
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
      fallbacks = sum(
        grepl("fallback", group$warning, ignore.case = TRUE)
      ),
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
      grid <- make_evaluation_grid(options$points)
      truth <- scenario$density(grid)
      finite_bandwidths <- log_bandwidth_grid(
        lower = 1e-4,
        upper = 0.25 - sqrt(.Machine$double.eps),
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
            kernel = "beta",
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
            kernel = "beta",
            start = "beta",
            bw = bandwidth
          )
        },
        scenario = scenario_name,
        replication = replication,
        estimator = "beta_start"
      )
    }
  }
}

results <- do.call(rbind, results)
summary <- summarize_results(results)
stem <- sprintf(
  "beta-mise-n%d-r%d-d%d-seed%d",
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
