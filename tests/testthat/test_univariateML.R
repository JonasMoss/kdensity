univariateML_densities <- c("norm", "lnorm", "weibull", "invgauss", "beta")

expect_true(density_namespace_available("dnorm"))
expect_true(density_namespace_available("stats::dnorm"))
expect_false(density_namespace_available("definitelymissingpkg::dfoo"))

methods::setClass("MockSupportKdensityCoverage", contains = "matrix", slots = c(type = "character"))
mock_support <- methods::new(
  "MockSupportKdensityCoverage",
  matrix(c(-Inf, Inf), nrow = 1),
  type = "R"
)
expect_equal(
  get_univariate_ml_support(list(support = mock_support)),
  list(bounds = c(-Inf, Inf), type = "R")
)

# There are some exceptions we won't test.

exceptions <- c("pareto", "unif", "lomax")

generate_random <- function(n, support) {
  if (support == "c(-Inf, Inf)") {
    return(rnorm(n))
  }
  if (support == "c(0, Inf)") {
    return(rgamma(n, 4, 4))
  }
  if (support == "c(1, Inf)") {
    return(rexp(n) + 1)
  }
  if (support == "c(0, 1)") {
    return(runif(n))
  }
  stop()
}

# Actual testing for all non-exceptions.
set.seed(10)
n <- 50

for (name in univariateML_densities) {
  support <- deparse(get_density_and_support(name)$support)
  rands <- generate_random(n, support)
  kde <- kdensity(rands, start = name)
  coef(kde)
  logLik(kde)
  AIC(kde)
}

# For Lomax
rands <- extraDistr::rlomax(n, 2, 3)
kde <- kdensity(rands, start = "lomax")
coef(kde)
logLik(kde)
AIC(kde)
