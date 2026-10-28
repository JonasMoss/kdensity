# Shared helpers for the gamma and beta MISE experiments. These functions
# evaluate fixed bandwidth grids across replications. The resulting minimizer
# of average ISE is an oracle-MISE bandwidth, not a sample-wise ISE oracle.

log_bandwidth_grid <- function(lower, upper, points,
                               include_infinite = FALSE) {
  if (points < 2L || lower <= 0 || upper <= lower) {
    stop("Oracle grids require points >= 2 and 0 < lower < upper.")
  }
  out <- exp(seq(log(lower), log(upper), length.out = points))
  if (include_infinite) {
    out <- c(out, Inf)
  }
  out
}

evaluate_oracle_grid <- function(x, evaluation_grid, truth, bandwidths,
                                 fit_function, scenario, replication,
                                 estimator) {
  rows <- vector("list", length(bandwidths))

  for (index in seq_along(bandwidths)) {
    bandwidth <- bandwidths[[index]]
    started <- proc.time()[["elapsed"]]
    error <- ""
    warning <- ""

    fit <- tryCatch(
      withCallingHandlers(
        fit_function(bandwidth),
        warning = function(condition) {
          warning <<- paste(
            c(warning, conditionMessage(condition))[nzchar(
              c(warning, conditionMessage(condition))
            )],
            collapse = " | "
          )
          invokeRestart("muffleWarning")
        }
      ),
      error = function(condition) {
        error <<- conditionMessage(condition)
        NULL
      }
    )

    ise <- NA_real_
    if (!is.null(fit)) {
      estimate <- tryCatch(
        fit(evaluation_grid),
        error = function(condition) {
          error <<- conditionMessage(condition)
          rep(NA_real_, length(evaluation_grid))
        }
      )
      if (all(is.finite(estimate))) {
        ise <- trapezoid_integral(
          evaluation_grid,
          (estimate - truth)^2
        )
      } else if (!nzchar(error)) {
        error <- "Density evaluation returned non-finite values."
      }
    }

    rows[[index]] <- data.frame(
      scenario = scenario,
      replication = replication,
      estimator = estimator,
      candidate = index,
      bandwidth = bandwidth,
      ise = ise,
      fit_seconds = proc.time()[["elapsed"]] - started,
      warning = warning,
      error = error,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

summarize_oracle_grid <- function(results) {
  groups <- split(
    results,
    interaction(
      results$scenario,
      results$estimator,
      results$candidate,
      drop = TRUE
    )
  )
  summaries <- lapply(groups, function(group) {
    successful <- is.finite(group$ise)
    ise <- group$ise[successful]
    data.frame(
      scenario = group$scenario[[1]],
      estimator = group$estimator[[1]],
      candidate = group$candidate[[1]],
      bandwidth = group$bandwidth[[1]],
      repetitions = nrow(group),
      successful = sum(successful),
      failures = sum(!successful),
      mise = if (length(ise)) mean(ise) else NA_real_,
      mise_se = if (length(ise) > 1L) {
        stats::sd(ise) / sqrt(length(ise))
      } else {
        NA_real_
      },
      mean_fit_seconds = mean(group$fit_seconds),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, summaries)
  summary[order(summary$scenario, summary$estimator, summary$candidate), ]
}

oracle_grid_optima <- function(summary) {
  groups <- split(
    summary,
    interaction(summary$scenario, summary$estimator, drop = TRUE)
  )
  optima <- lapply(groups, function(group) {
    valid <- is.finite(group$mise)
    if (!any(valid)) {
      return(group[1L, , drop = FALSE])
    }
    group[which.min(replace(group$mise, !valid, Inf)), , drop = FALSE]
  })
  out <- do.call(rbind, optima)
  rownames(out) <- NULL
  out[order(out$scenario, out$estimator), ]
}

global_series_test_statistic <- function(selector, n) {
  indices <- 3:length(selector$theta)
  coefficients <- unname(selector$theta[indices])
  covariance <- selector$coefficient_covariance[
    indices,
    indices,
    drop = FALSE
  ]
  decomposition <- eigen(covariance, symmetric = TRUE)
  keep <- decomposition$values >
    max(decomposition$values) * sqrt(.Machine$double.eps)

  if (!any(keep)) {
    statistic <- 0
    degrees_of_freedom <- 0L
    p_value <- 1
  } else {
    coordinates <- drop(
      crossprod(decomposition$vectors[, keep, drop = FALSE], coefficients)
    )
    statistic <- n * sum(coordinates^2 / decomposition$values[keep])
    degrees_of_freedom <- sum(keep)
    p_value <- stats::pchisq(
      statistic,
      df = degrees_of_freedom,
      lower.tail = FALSE
    )
  }

  list(
    statistic = statistic,
    degrees_of_freedom = degrees_of_freedom,
    p_value = p_value
  )
}

conservative_series_alpha <- function(n) {
  # A fixed alpha retains a fixed probability of perturbing an exact
  # parametric model. This sequence removes that asymptotically while still
  # allowing consistent tests to detect fixed alternatives.
  min(0.01, 1 / n)
}

apply_global_series_test <- function(selector, n, alpha = 0.05) {
  test <- global_series_test_statistic(selector, n)
  selector$unregularized_bandwidth <- selector$bandwidth
  selector$global_statistic <- test$statistic
  selector$global_degrees_of_freedom <- test$degrees_of_freedom
  selector$global_p_value <- test$p_value
  selector$global_asymptotic_p_value <- test$p_value
  selector$global_selected <- test$p_value < alpha
  if (!selector$global_selected) {
    selector$bandwidth <- Inf
  }
  selector
}

apply_bootstrap_global_series_test <- function(
    selector,
    n,
    simulate,
    refit,
    replicates = 199L,
    alpha = conservative_series_alpha(n)) {
  replicates <- as.integer(replicates)
  if (length(replicates) != 1L || is.na(replicates) || replicates < 19L) {
    stop("'replicates' must be a single integer of at least 19.")
  }
  if (!is.function(simulate) || !is.function(refit)) {
    stop("'simulate' and 'refit' must be functions.")
  }

  observed <- global_series_test_statistic(selector, n)
  selector$unregularized_bandwidth <- selector$bandwidth
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    saved_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", saved_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  bootstrap_statistics <- rep(NA_real_, replicates)
  for (replication in seq_len(replicates)) {
    bootstrap_statistics[[replication]] <- tryCatch(
      {
        bootstrap_x <- simulate(n)
        bootstrap_selector <- refit(bootstrap_x)
        global_series_test_statistic(
          bootstrap_selector,
          length(bootstrap_x)
        )$statistic
      },
      error = function(condition) NA_real_
    )
  }

  successful <- is.finite(bootstrap_statistics)
  successful_replicates <- sum(successful)
  if (successful_replicates < ceiling(0.8 * replicates)) {
    stop(
      "Too few successful parametric-bootstrap refits: ",
      successful_replicates,
      " of ",
      replicates,
      "."
    )
  }
  exceedances <- sum(
    bootstrap_statistics[successful] >= observed$statistic
  )
  p_value <- (1 + exceedances) / (1 + successful_replicates)
  critical_value <- unname(stats::quantile(
    bootstrap_statistics[successful],
    probs = 1 - alpha,
    names = FALSE,
    type = 1
  ))

  selector$global_statistic <- observed$statistic
  selector$global_degrees_of_freedom <- observed$degrees_of_freedom
  selector$global_asymptotic_p_value <- observed$p_value
  selector$global_p_value <- p_value
  selector$global_selected <- p_value <= alpha
  selector$bootstrap_replicates <- replicates
  selector$bootstrap_successful <- successful_replicates
  selector$bootstrap_failures <- replicates - successful_replicates
  selector$bootstrap_exceedances <- exceedances
  selector$bootstrap_alpha <- alpha
  selector$bootstrap_critical_value <- critical_value
  if (!selector$global_selected) {
    selector$bandwidth <- Inf
  }
  selector
}

bootstrap_transition_weight <- function(selector, power = 4) {
  if (length(power) != 1L || !is.finite(power) || power <= 0) {
    stop("'power' must be a single positive number.")
  }
  statistic <- selector$global_statistic
  critical_value <- selector$bootstrap_critical_value
  if (!is.finite(statistic) || statistic <= 0) {
    return(0)
  }
  if (!is.finite(critical_value) || critical_value <= 0) {
    return(1)
  }
  stats::plogis(power * (log(statistic) - log(critical_value)))
}

smooth_bootstrap_bandwidth <- function(selector, power = 4) {
  weight <- bootstrap_transition_weight(selector, power)
  bandwidth <- selector$unregularized_bandwidth
  selector$smooth_weight <- weight
  selector$transition_power <- power
  selector$effective_roughness <- weight * selector$roughness
  selector$bandwidth <- if (!is.finite(bandwidth) || weight == 0) {
    Inf
  } else {
    bandwidth / weight^(2 / 5)
  }
  selector
}

weighted_density_fit <- function(parametric_fit, nonparametric_fit, weight) {
  if (length(weight) != 1L || !is.finite(weight) ||
      weight < 0 || weight > 1) {
    stop("'weight' must be a single number in [0, 1].")
  }
  bw <- nonparametric_fit$bw
  out <- function(y) {
    (1 - weight) * parametric_fit(y) + weight * nonparametric_fit(y)
  }
  class(out) <- "kdensity"
  out
}

select_series_degree <- function(selector, n, penalty = log(n)) {
  maximum_degree <- selector$degree
  scores <- numeric(maximum_degree + 1L)
  names(scores) <- 0:maximum_degree
  scores[1L] <- 0
  scores[2L] <- -Inf

  for (order in 2:maximum_degree) {
    indices <- 3:(order + 1L)
    coefficients <- unname(selector$theta[indices])
    covariance <- selector$coefficient_covariance[
      indices,
      indices,
      drop = FALSE
    ]
    decomposition <- eigen(covariance, symmetric = TRUE)
    keep <- decomposition$values >
      max(decomposition$values) * sqrt(.Machine$double.eps)
    coordinates <- drop(
      crossprod(decomposition$vectors[, keep, drop = FALSE], coefficients)
    )
    statistic <- n * sum(coordinates^2 / decomposition$values[keep])
    scores[[order + 1L]] <- statistic - sum(keep) * penalty
  }

  selected_degree <- which.max(scores) - 1L
  list(
    degree = selected_degree,
    scores = scores
  )
}
