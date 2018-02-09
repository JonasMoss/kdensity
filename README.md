kdensity: Kernel density estimation with a parametric start
================
Jonas Moss
8 February 2018

## Introduction

Kernel density estimation with a parametric start was introduced by Nils
Lid Hjort and Ingrid Glad in their 1995 paper [Nonparametric Density
Estimation with a Parametric
Start](https://projecteuclid.org/euclid.aos/1176324627). The idea is to
start out with a parametric density before you do your kernel density
estimation, so that your actual kernel density estimation will be a
correction to the original parametric estimate. Why is this a good idea?
Because the resulting estimator will be better than a naïve kernel
density estimator whenever the true density is close to your suggestion,
and probably not worse if not. Read the paper for more information.

The goal of this `R` package is to make kernel density estimation with a
parametric start both routine and flexible. In order to achieve this, we
have decided on emulating the behaviour of the function `density` in the
`stats` package (which is included in your `R` installation) to the
highest extent possible.

In addition to parametric starts, the package implements some
*asymmetric kernels*. These kernels are useful when modelling data with
sharp boundaries, such as data supported on the positive half-line or
the unit interval. Currently we support the following asymmetric
kernels:

  - MC Jones and Daniel Henderson’s Gaussian copula KDE, from
    [Kernel-Type Density Estimation on the Unit Interval
    (2007)](https://academic.oup.com/biomet/article-abstract/94/4/977/246269).
    This is used for data on the unit interval. The bandwidth selection
    mechanism described in that paper is implemented as well. This
    kernel is called `gcopula`.

  - Song Xi Chen’s two gammakernels from [Probability Density Function
    Estimation Using Gamma Kernels
    (2000)](https://link.springer.com/article/10.1023/A:1004165218295).
    This is used for data supported on the positive half-line. These
    kernels are called `gamma` and `gamma_biased`.

So the package has two main features: parametric starts and asymmetric
kernels. The features can be combined to make asymmetric kernel
densities estimators with parametric starts, see the example below. The
package contains only one function, `kdensity`, in addition to the
generics `plot`, `points`, `lines`, `summary`, and `print`. It’s not
completely documentet yet, so stay tuned for that.

## Example

The data set `sunspot.month` is described at `R` site
[here](https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/sunspot.month.html).
This is a good example to showcase the usefulness of kernel density
estimation with a parametric start, as it is supported on positive
half-line. In such a case, traditional density estimators are prone to
serious *boundary bias*. However, parametric starts allows us to
circumvent this by using a pre-specified parametric density supported on
the half-line, for instance the exponential, which clearly fits the data
well.

``` r
library("kdensity")
data = sunspot.month[sunspot.month != 0]
x = seq(0, 300, by = 0.5)
hist(data, freq = FALSE, breaks = 40,
     main = "Monthly Sunspot Numbers, 1749 – 1983",
     xlab = "Monthly mean relative sunspot numbers")
lines(kdensity(data, start = "exponential", kernel = "gaussian", adjust = 1),
      lwd = 2, lty = 2, col = "red")
lines(x, dexp(x, 1/mean(data)), lwd = 2, lty = 3, col = "blue")
```

<img src="README_files/figure-gfm/gaussiankernelexp plot-1.png" width="750px" />

As seen from the histogram, the fit is quite bad, so let us try out the
`gamma` kernel. Notice the argument `start = "uniform"`, which states
that the parametric start is *constant*. Using this argument reduces the
density estimates to a classical kernel density estimate.

``` r
library("kdensity")
data = sunspot.month[sunspot.month != 0]
x = seq(0, 300, by = 0.5)
hist(data, freq = FALSE, breaks = 40,
     main = "Monthly Sunspot Numbers, 1749 – 1983",
     xlab = "Monthly mean relative sunspot numbers")
lines(kdensity(data, start = "uniform", kernel = "gamma", adjust = 1),
      lwd = 2, lty = 2, col = "red")
lines(x, dexp(x, 1/mean(data)), lwd = 2, lty = 3, col = "blue")
```

<img src="README_files/figure-gfm/gammakernel plot-1.png" width="750px" />

While this helped, the fit could certainly be better, as the density
estimate is artificially pushed downwards at the boundary. So let us try
a `gamma` kernel with an `exponential` start:

``` r
library("kdensity")
data = sunspot.month[sunspot.month != 0]
x = seq(0, 300, by = 0.5)
hist(data, freq = FALSE, breaks = 40,
     main = "Monthly Sunspot Numbers, 1749 – 1983",
     xlab = "Monthly mean relative sunspot numbers")
lines(kdensity(data, start = "exponential", kernel = "gamma", adjust = 1),
      lwd = 2, lty = 2, col = "red")
lines(x, dexp(x, 1/mean(data)), lwd = 2, lty = 3, col = "blue")
```

<img src="README_files/figure-gfm/gammakernelexp plot-1.png" width="750px" />

This fit is much better now. What’s more, it gives qualitatively
different predictions about the probability of extremely many sunspots
than the exponential does.

## How to use

The function `kdensity` takes some `data`, a kernel `kernel` and a
parametric start `start`. Currently the following arguments for a
symmetric `kernel` are supported: `gaussian`, `laplace`, `epanechnikov`,
`rectangular`, `triangular`, `biweight`, `triweight`, `tricube`,
`cosine`, `optcosine`. In addition, `gcopula` is supported for the unit
interval, while `gamma` and `gamma_biased` are supported for the positve
half-line. For the `start` parameter, the following are accepted (as
strings\!): `uniform`, `normal`, `gamma`, `exponential`,
`inverse_gaussian`, `lognormal`, `beta`, and `laplace`. Their parameters
are estimated by maximum likelihood. In addition to this, you can
specify the `support` parameter, which is used to find the normalizing
constant.You do not have to worry about this one though, as it is
automatically deduced from the choice of `start`. Still, if you want to
have support on \[0, 12.3\], it is possible.

Finally, instead of specifying `start` as a string, you can send a
`list` containing three named elements: `density` is a density function
with named parameters (where `x` is the main argument), `estimator` is a
function of `data` that estimates the parameters for you, and `support`
is the domain of defition for the density. For example, take a look at
this implementation of the log-normal distribution:

``` r
lognormal = list(
  density = dlnorm,
  estimator = function(data) {
    c(meanlog = mean(log(data)),
      sdlog = sd(log(data)))
  },
  support = c(0, Inf)
)
```

Using this definition we can run a new `kdensity` function as follows:

``` r
plot(kdensity(data, start = lognormal, kernel = "gamma"), 
     lwd = 2, main = "Log normal plot")
```

<img src="README_files/figure-gfm/lnormplot-1.png" width="750px" />

The `plot` function works just as in the case of `stats::density`.
Moreover, `lines` and `points` does as well. Since the return value of
`kdensity` is a function, it is callable, as in:

``` r
new_density = kdensity(data, start = "exponential")
new_density(56.7)
```

    ## [1] 0.008375948

## Installation

First you need to install the package `devtools` from `CRAN`. From
inside `R`, use the following command.

``` r
devtools::install_github("JonasMoss/kdensity")
```

This installs the latest version of the package from GitHub. Enjoy, and
please feel compelled to e-mail us suggestions and bug reports.
