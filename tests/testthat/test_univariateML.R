univariateML_densities <- intersect(
  names(starts),
  univariateML::univariateML_models
)

univariateML_densities <- c("norm", "lnorm", "weibull", "invgauss", "beta")

# There are some exceptions we won't test.

exceptions <- c("pareto", "unif", "lomax")

univariateML_densities <- setdiff(
  univariateML_densities,
  exceptions
)

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
  univariateML::univariateML_metadata[[paste0("ml", name)]]
}

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
