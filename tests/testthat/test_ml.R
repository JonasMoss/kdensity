context("maximum likelihood")

## The test MLEs are obtained from nlm.

### Checking beta:

set.seed(313)
data1 = runif(100)
data2 = rbeta(10, 2, 7)

mle1 = nlm(function(p) {
  -mean(dbeta(data1, p[1], p[2], log = TRUE))
  }, p = c(1,1))$estimate

mle2 = nlm(function(p) {
  -mean(dbeta(data2, p[1], p[2], log = TRUE))
}, p = c(1,1))$estimate

names(mle1) = c("shape1", "shape2")
names(mle2) = c("shape1", "shape2")

expect_equal(mle1, mlbeta(data1), tolerance = 1e-5)
expect_equal(mle2, mlbeta(data2), tolerance = 1e-5)
expect_equal(mlbeta(data2, type = "gradient"), mlbeta(data2), tolerance = 1e-5)
expect_equal(mlbeta(data1, type = "gradient"), mlbeta(data1, type = "hessian"), tolerance = 1e-5)


### Checking gamma:

set.seed(313)
data1 = rgamma(100, 1, 1)
data2 = rgamma(10, 3, 7)

mle1 = nlm(function(p) {
  -mean(dgamma(data1, p[1], p[2], log = TRUE))
}, p = c(1,1))$estimate

mle2 = nlm(function(p) {
  -mean(dgamma(data2, p[1], p[2], log = TRUE))
}, p = c(1,1))$estimate

names(mle1) = c("shape", "rate")
names(mle2) = c("shape", "rate")

expect_equal(mle1, mlgamma(data1), tolerance = 1e-5)
expect_equal(mle2, mlgamma(data2), tolerance = 1e-5)
expect_warning(mlgamma(data2, iterlim = 1))

### Checking Weibull:

set.seed(313)
data1 = rweibull(100, 1, 1)
data2 = rweibull(10, 3, 7)

mle1 = nlm(function(p) {
  -mean(dweibull(data1, p[1], p[2], log = TRUE))
}, p = c(1,1))$estimate

mle2 = nlm(function(p) {
  -mean(dweibull(data2, p[1], p[2], log = TRUE))
}, p = c(3,7))$estimate

names(mle1) = c("shape", "scale")
names(mle2) = c("shape", "scale")

expect_equal(mle1, mlweibull(data1), tolerance = 1e-5)
expect_equal(mle2, mlweibull(data2), tolerance = 1e-5)
expect_warning(mlweibull(data2, iterlim = 1))

### Checking gumbel:

set.seed(313)
data1 = extraDistr::rgumbel(100, 1, 1)
data2 = extraDistr::rgumbel(10, 3, 7)

mle1 = nlm(function(p) {
  -mean(extraDistr::dgumbel(data1, p[1], p[2], log = TRUE))
}, p = c(1,1))$estimate

mle2 = nlm(function(p) {
  -mean(extraDistr::dgumbel(data2, p[1], p[2], log = TRUE))
}, p = c(3,7))$estimate

names(mle1) = c("loc", "scale")
names(mle2) = c("loc", "scale")

expect_equal(mle1, mlgumbel(data1), tolerance = 1e-5)
expect_equal(mle2, mlgumbel(data2), tolerance = 1e-5)
expect_warning(mlgumbel(data2, iterlim = 1))


### Checking Kumaraswamy:

set.seed(313)
data1 = extraDistr::rkumar(100, 4, 4)
data2 = extraDistr::rkumar(10, 3, 7)

mle1 = nlm(function(p) {
  -mean(extraDistr::dkumar(data1, p[1], p[2], log = TRUE))
}, p = c(4, 4))$estimate

mle2 = nlm(function(p) {
  -mean(extraDistr::dkumar(data2, p[1], p[2], log = TRUE))
}, p = c(3, 7))$estimate

names(mle1) = c("a", "b")
names(mle2) = c("a", "b")

expect_equal(mle1, mlkumar(data1), tolerance = 1e-4)
expect_equal(mle2, mlkumar(data2), tolerance = 1e-4)
expect_warning(mlkumar(data2, iterlim = 1))


