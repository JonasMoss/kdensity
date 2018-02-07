kdensity: Kernel density estimation with a parametric start
================
Jonas Moss
7 February 2018

## Introduction

Kernel density estimation with a parametric start was introduced by Nils
Lid Hjort and Ingrid Glad in their 1996 paper [Nonparametric Density
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

**NB:** This package is mint out of box\! At least one core feature, the
choice of bandwidth, has not been implemented, and the documentation is
not fully spelled out. Come back later for more\! The package only
supports estimates in one dimension at the moment.

## Example

The data set `sunspot.month` is described at `R` site
[here](https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/sunspot.month.html).
This is a good example to showcase the usefulness of kernel density
estimation with a parametric start, as it is supported on positive
half-line. In such a case, traditional density estimators are prone to
serious *boundary bias*. However, parametric starts allows us to
circumvent this by using a pre-specified parametric density supported on
the half-line, for instance the
exponential.

``` r
plot(kdensity::kdensity(sunspot.month, start = "exponential", bw = 16), lwd = 2)
```

<img src="README_files/figure-gfm/small plot-1.png" width="750px" />

Lets compare this to the behaviour of `density` and `dexp`:

``` r
y = seq(0, max(sunspot.month) + 10, by = 0.5)
plot(kdensity::kdensity(sunspot.month, start = "exponential", bw = 16), lwd = 2)
lines(y, dexp(y, 1/mean(sunspot.month)), lty = 2)
lines(density(sunspot.month, kernel = "gauss", adjust = 2), lty = 3)
```

<img src="README_files/figure-gfm/full plot-1.png" width="750px" /> Here
`density` is easily seen to be far more affected by boundary bias than
`kdensity`, and it capturs more of features of the data than the
exponential distribution does.

## How to use

The function `kdensity` takes some `data`, a kernel `kernel` and a
parametric start `start`. Currently the following arguments for `kernel`
are supported: `gaussian`, `laplace`, `epanechnikov`, `rectangular`,
`triangular`, `biweight`, `triweight`, `tricube`, `cosine`, and
`optcosine`. For the `start` parameter, the following are accepted (as
strings\!): `uniform`, `normal`, `gamma`, `exponential`,
`inverse_gaussian`, `lognormal`, `beta`, and `laplace`. Their parameters
are estimated by maximum likelihood. In addition to this, you can
specify the `support` parameter, which is used to find the normalizing
constant.You do not have to worry about this one though, as it is
automatically deduced from the choice of `start`. Still, if you want to
have support on \[0, 12.3\], it is possible.

Finally, instead of specifying `start` as a string, you can send a
`list` containing two named functions as elements: `density` is a
density function with named parameters, while `estimator` is a function
of `data` that estimates the parameters for you. For example, take a
look at this implementation of the inverse Gaussian distribution:

``` r
start_inverse_gaussian = list(
  density = statmod::dinvgauss,
  estimator = function(data) {
    c(mean       = mean(data),
      dispersion = mean(1/data - 1/mean(data)))
  }
)
```

Using this definition we can run a new `kdensity` function as
follows:

``` r
plot(kdensity::kdensity(sunspot.month + 1, start = start_inverse_gaussian, support = c(0, Inf), bw = 16), lwd = 2)
```

<img src="README_files/figure-gfm/lnormplot-1.png" width="750px" /> Here
`sunspot.month + 1` is used inside `kdensity` in order to make the data
strictly positive.

The `plot` function works just as in the case of `stats::density`.
Moreover, `lines` and `points` does as well. Since the return value of
`kdensity` is a function, it is callable, as
in:

``` r
new_density = kdensity::kdensity(sunspot.month + 1, start = "lognormal", bw = 16)
new_density(56.7)
```

    ## [1] 0.007617428

## Installation

First you need to install the package `devtools` from `CRAN`. Then you
can call:

``` r
devtools::install_github("JonasMoss/kdensity")
```

Enjoy\!
