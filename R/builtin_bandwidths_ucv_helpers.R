kdensity_sq <- function(x, h, kernel_fun, parametric_start, parametric_start_data,
                        parametric_start_vector, support) {
  normalization <- 1

  pre_function <- function(y) {
    sapply(y, function(y) mean(1 / h * kernel_fun(y, x, h) / parametric_start_data) * parametric_start_vector(y))
  }

  normalization <- tryCatch(stats::integrate(pre_function, lower = support[1], upper = support[2])$value,
    error = function(e) {
      stop("Normalization error: The function will not integrate. Two common causes are: 1.) The kernel is non-smooth, try a smooth kernel if possible. 2.) The supplied support is incorrect.")
    }
  )

  return_function <- function(y) {
    n <- length(y)
    parametric_start_vector_y <- parametric_start_vector(y)
    sapply(1:n, function(i) {
      (1 / h * mean(kernel_fun(y[i], x, h) * parametric_start_vector_y[i] / parametric_start_data) / normalization)^2
    })
  }
  return_function
}
