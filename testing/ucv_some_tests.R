library("magrittr")

## Gumbel and mtcars$mpg
data = mtcars$mpg
kdensity(data, kernel = "normal", start = "gumbel", bw = "ucv") %>%
  plot(main = "Miles per Gallon") %>%
  lines(col = "red", plot_start = TRUE)
rug(data)


## Weibull and ucv.
set.seed(313)
data = rweibull(100, 2, 7)
kdensity(data, kernel = "gaussian", start = "normal", bw = "ucv", support = c(0, Inf)) %>%
  plot %>%
  lines(col = "red", plot_start = TRUE)
rug(data)

## Gamma kernels and ucv.
xx = seq(0, 20, by = 0.01)
data = rbeta(100, 2, 7)
kdensity(data, kernel = "gamma", start = "uniform", bw = "ucv") %>%
  plot(col = "blue")
kdensity(data, kernel = "gamma", start = "gumbel", bw = "ucv") %>%
  lines(col = "red")
rug(data)
lines(xx, dweibull(xx, 2, 7), col = "purple")

## Betas and ucv.
set.seed(313)
xx = seq(0, 1, by = 0.01)
kdensity(data, kernel = "beta", start = "uniform", bw = "ucv") %>%
  plot(col = "blue")
kdensity(data, kernel = "beta", start = "beta", bw = "ucv") %>%
  lines(col = "red") %>%
  lines(col = "black", plot_start = TRUE)
rug(data)
lines(xx, dgamma(xx, 2, 7), col = "purple")

## Copulas and ucv.
kdensity(data, kernel = "gcopula", start = "uniform", bw = "ucv") %>%
  plot(col = "blue")
kdensity(data, kernel = "gcopula", start = "beta", bw = "ucv") %>%
  lines(col = "red") %>%
  lines(col = "black", plot_start = TRUE)
rug(data)
lines(xx, dgamma(xx, 2, 7), col = "purple")


## Pareto and ucv.
set.seed(313)
data = 1/(runif(100))^(1/3)
xx = seq(0, 20, by = 0.01)
kdensity(data, kernel = "gamma", start = "pareto", bw = "ucv") %>%
  plot(col = "blue")
kdensity(data, kernel = "normal", start = "pareto", bw = "ucv") %>%
  lines(col = "red") %>%
  lines(col = "black", plot_start = TRUE)
rug(data)

