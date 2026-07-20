#' @useDynLib kdensity, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

make_kernel_evaluator <- function(x, h, kernel, parametric_start_data) {
  weights <- 1 / (length(x) * h * parametric_start_data)
  kernel_id <- kernel$compiled_id

  if (is.null(kernel_id)) {
    return(function(y) {
      vapply(y, function(point) {
        sum(kernel$kernel(point, x, h) * weights)
      }, numeric(1))
    })
  }

  cache1 <- numeric()
  cache2 <- numeric()
  sorted_supported <- !is.null(kernel_id) && kernel_id %in% 2:10 &&
    is.finite(h) && h > 0 &&
    all(is.finite(x)) && all(is.finite(weights))
  sorted_ready <- FALSE
  sorted_x <- numeric()
  sorted_weights <- numeric()
  sorted_cache <- list()
  cumulative_work <- 0
  # Building the compact moment tree costs about as much as 150 direct
  # evaluations per observation. Laplace only needs sorting and pays back
  # after roughly ten. Accumulate work across scalar calls as well as batches.
  sorted_work_threshold <- if (identical(kernel_id, 10L)) {
    max(10000, 10 * length(x))
  } else {
    max(50000, 150 * length(x))
  }

  if (kernel_id == 11L) {
    cache1 <- stats::qnorm(x)
  } else if (kernel_id %in% 12:15) {
    cache1 <- log(x)
    if (kernel_id %in% 14:15) {
      cache2 <- log1p(-x)
    }
  }

  function(y) {
    if (sorted_supported) {
      cumulative_work <<- cumulative_work +
        as.double(length(x)) * as.double(length(y))
      if (sorted_ready || cumulative_work >= sorted_work_threshold) {
        if (!sorted_ready) {
          indices <- order(x, method = "radix")
          sorted_x <<- x[indices]
          sorted_weights <<- weights[indices]
          if (kernel_id != 10L) {
            sorted_cache <<- compact_kernel_cache_cpp(
              sorted_x,
              sorted_weights,
              h
            )
          }
          sorted_ready <<- TRUE
        }
        return(kernel_weighted_sum_sorted_cpp(
          y = y,
          sorted_x = sorted_x,
          sorted_weights = sorted_weights,
          h = h,
          kernel_id = kernel_id,
          cache = sorted_cache
        ))
      }
    }

    kernel_weighted_sum_cpp(
      y = y,
      x = x,
      weights = weights,
      h = h,
      kernel_id = kernel_id,
      cache1 = cache1,
      cache2 = cache2
    )
  }
}

compiled_start_parameter_names <- function(start_id) {
  switch(as.character(start_id),
    "1" = c("mean", "sd"),
    "2" = "rate",
    "3" = c("meanlog", "sdlog"),
    "4" = c("shape", "rate"),
    "5" = c("shape", "scale"),
    "6" = c("shape1", "shape2"),
    NULL
  )
}

make_compiled_loo_parameters <- function(x, start) {
  start_id <- start$compiled_id
  parameter_names <- compiled_start_parameter_names(start_id)
  if (is.null(start_id) || is.null(parameter_names)) {
    return(NULL)
  }

  parameters <- designated_loo_parameters_cpp(x, start_id)
  colnames(parameters) <- parameter_names
  parameters
}

compiled_ucv_term2 <- function(x, h, kernel, start, parameters) {
  kernel_id <- kernel$compiled_id
  start_id <- start$compiled_id
  if (is.null(kernel_id) || is.null(start_id) || is.null(parameters)) {
    return(NULL)
  }

  ucv_loo_term_cpp(
    x = x,
    h = h,
    kernel_id = kernel_id,
    start_id = start_id,
    parameters = parameters
  )
}

known_normalization <- function(start, kernel, support, h, x, estimates) {
  start_id <- start$normalization_id
  kernel_id <- kernel$compiled_id
  if (is.null(start_id) || is.null(kernel_id)) {
    return(NULL)
  }

  if (start_id == 1L && kernel_id %in% 1:10 && h > 0 &&
      identical(support, c(-Inf, Inf))) {
    return(1)
  }

  if (start_id == 1L && kernel_id == 11L &&
      identical(support, c(0, 1)) &&
      h > 0 && h < sqrt(2) &&
      all(x > 0 & x < 1)) {
    rho <- 1 - h^2
    return(sqrt(1 - rho^2) / h)
  }

  if (start_id == 2L && kernel_id == 1L &&
      identical(support, c(-Inf, Inf))) {
    if (!all(c("mean", "sd") %in% names(estimates))) {
      return(NULL)
    }

    mean_value <- unname(estimates[["mean"]])
    sd_value <- unname(estimates[["sd"]])
    if (!is.finite(mean_value) || !is.finite(sd_value) ||
        sd_value <= 0 || h <= 0) {
      return(NULL)
    }

    log_ratio <- stats::dnorm(
      x,
      mean = mean_value,
      sd = sqrt(sd_value^2 + h^2),
      log = TRUE
    ) - stats::dnorm(
      x,
      mean = mean_value,
      sd = sd_value,
      log = TRUE
    )
    return(mean(exp(log_ratio)))
  }

  NULL
}
