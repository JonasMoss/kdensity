library("magrittr")

## Testing custom start:
skew_hyperbolic = list(
  density   = SkewHyperbolic::dskewhyp,
  estimator = function(x) SkewHyperbolic::skewhypFit(x, printOut = FALSE)$param,
  support   = c(-Inf, Inf)
)

# Gumbel and mtcars$mpg
data = mtcars$mpg
kdensity(data, kernel = "normal", start = skew_hyperbolic, bw = "ucv") %>%
  plot(main = "Miles per Gallon") %>%
  lines(col = "red", plot_start = TRUE) %>%
  coef
rug(data)

## Testing custom kernel:

t4 = list(
  kernel   = function(y, x, h) dt((x-y)/h, df = 4),
  sd       = 1.41,
  support   = c(-Inf, Inf)
)


# Gumbel and mtcars$mpg
data = mtcars$mpg
kdensity(data, kernel = t4, start = "gumbel", bw = "nrd0") %>%
  plot(main = "Miles per Gallon") %>%
  lines(col = "red", plot_start = TRUE) %>%
  coef
rug(data)
