---
title: 'kdensity: An R package for kernel density estimation with parametric starts and asymmetric Kernels'
tags:
  - R
  - statistics
  - kernel density estimation
  - non-parametric statistics
  - non-parametrics
  - non-parametric density estimation
  - boundary bias
authors:
  - name: Jonas Moss
    orcid: 0000-0002-6876-6964
    affiliation: 1
  - name: Martin Tveten
    orcid: 0000-0002-4236-633X
    affiliation: 1
affiliations:
 - name: University of Oslo
   index: 1
date: 11 July 2019
bibliography: paper.bib
---

# Summary

Kernel density estimation [@silverman2018density] is a popular method for 
non-parametric density estimation based on placing kernels on each data point. 
@hjort_glad_1995 extended kernel density estimation with *parametric starts*.
The parametric start is a parametric density that is multiplied with the kernel
estimate. When the data-generating density is reasonably close to the parametric
start density, kernel density estimation with that parametric start will outperform
ordinary kernel density estimation.

Asymmetric kernels are useful for estimating densities on the half-open interval $\left[0,\infty\right)$ and bounded intervals such as $\left[0, 1\right]$. On 
such intervals symmetric kernels are prone to serious boundary bias that should
be corrected [@marron1994transformations]. Asymmetric kernels are designed to
avoid boundary bias.

`kdensity` is an R package [@r] to calculate and display kernel density 
estimates using non-parametric starts and potentially asymmetric kernels. In 
addition to the classical symmetric kernels, `kdensity` supports the following 
asymmetric kernels: For the unit interval, the Gaussian copula kernel of @jones2007miscellanea and the beta kernels of @chen1999beta are supported. On 
the half-open interval the gamma kernel of @chen2000probability is supported. 
The supported non-parametric starts include the normal, Laplace, Gumbel, 
exponential, gamma, log-normal, inverse Gaussian, Weibull, Beta, and Kumaraswamy
densities. The parameters of all parametric starts are estimated using maximum 
likelihood. The implemented bandwidth selectors are the classical bandwidth 
selectors from `stats`, unbiased cross-validation, the Hermite polynomial method 
from @hjort_glad_1995, and the tailored bandwidth selector for the Gaussian 
copula method of @jones2007miscellanea. User defined parametric starts, 
kernels and bandwidth selectors are also supported. 

The following example uses the \code{airquality} data set from the built-in
R package `datasets`. Since the data is positive we use Chen's gamma kernel. 
As the data is likely to be better approximated by a gamma distribution than a 
uniform distribution, we use the gamma parametric start. The plotted density is
in figure 1, where the gamma distribution with parameters estimated by maximum 
likelihood is in red and the ordinary kernel density estimate in blue. 
Notice the boundary bias of the ordinary kernel density estimator. 

```r
# install.packages("kdensity")
library("kdensity")
kde = kdensity(airquality$Wind, start = "gamma", kernel = "gamma")
plot(kde, main = "Wind speed (mph)")
lines(kde, plot_start = TRUE, col = "red")
lines(density(airquality$Wind, adjust = 2), col = "blue")
rug(airquality$Wind)
```
![The *airquality* data set. Kernel density estimate in black and estimated gamma distribution in red.](example.png)

# References
