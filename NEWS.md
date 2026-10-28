## kdensity 2.0.0 (2026-10-28)

kdensity is now in maintenance mode: future releases fix bugs and keep the
package working on CRAN, but add no new features.

### Breaking changes

* `kernel = "gamma"` now implements Chen's (2000) modified gamma kernel, a
  gamma density with scale `h` and shape `y / h` for `y >= 2h` and
  `(y / (2h))^2 + 1` below. Earlier versions used a shape one larger, so for
  `y >= 2h` the kernel coincided with the uncorrected `gamma_biased` kernel and
  its boundary modification increased the error. Every fit that uses the gamma
  kernel changes, including the bandwidths selected for it.

* Every built-in kernel now defaults to a closed-form bandwidth rule, so no
  default needs cross-validation. A gamma-kernel fit at n = 10,000 took
  48--198 seconds with the old `"ucv"` default; it now takes milliseconds.
  The new defaults:
  - `gamma` and `gamma_biased`: `"GR"` (was `"ucv"`).
  - `beta` and `beta_biased`: `"HS"` (was `"ucv"`, except for `beta` with a
    uniform start).
  - `gcopula`: `"SJ"` (was `"JH"` with a uniform start and `"ucv"`
    otherwise). `"JH"` oversmooths multimodal data badly, up to 29 times the
    optimal integrated squared error in simulations.

  The rules ignore the parametric start, as `"RHE"` already did for symmetric
  kernels, so they undersmooth a little when the start fits well; `bw = "ucv"`
  remains available for a start-aware choice. In simulations over gamma,
  exponential, lognormal, Weibull, beta and mixture targets, the new defaults
  did as well as or better than `"ucv"` with uniform starts.

* `bw = "SJ"` with the `gcopula` kernel now applies Sheather and Jones's
  selector on the probit scale, where the Gaussian copula estimator is a
  Gaussian kernel estimate, and converts the result to the copula bandwidth.
  It previously applied the selector to the raw data, which does not give a
  bandwidth for this kernel.

* Fitted `kdensity` objects are read-only. The `update()`, `$<-` and `[[<-`
  methods are removed, so `kde$bw <- 0.5` now gives an error. `update()` falls
  back to `stats::update.default()`, which refits from the stored call with
  the changed arguments (as for `lm()`) and returns a new object: write
  `kde <- update(kde, bw = 0.5)`. This also makes `update()` work for fits with
  a numeric or custom-function bandwidth. The data must be visible from where
  `update()` is called.

* Custom kernels and parametric starts must use names that differ from the
  built-ins (for example `my_gaussian`). Before, a custom kernel or start with a
  built-in name silently replaced it for the rest of the session and disabled
  its compiled path. Re-registering a custom name is still allowed.

### New features

* `kdensity_grid()` evaluates a fit on a grid by linear binning and FFT
  convolution, for every built-in symmetric kernel and for the Gaussian-copula
  kernel (on the probit scale). It convolves with the fitted kernel itself, so
  it also covers kernels that `stats::density()` lacks (triweight, tricube,
  laplace). The result is approximate; calling a `kdensity` object stays exact.

* `plot()`, `lines()` and `points()` gain a `grid` argument. The default,
  `"auto"`, uses FFT for the Gaussian and Gaussian-copula kernels and exact
  evaluation for the others, where it is as fast. Use `grid = "exact"` or
  `grid = "fft"` to choose.

* New bandwidth selector `"GR"`: Chen's AMISE-optimal bandwidth for the
  bias-corrected gamma kernel, with a fitted gamma density in place of the
  unknown one. The fitted shape is floored at 3 (keeping the fitted mean),
  since the rule needs a shape above 3/2, and exact zeros are left out of the
  reference fit.

### Bug fixes

* `bw = "ucv"` for asymmetric kernels and for parametric starts compared the
  integral of the squared *normalized* estimate with *unnormalized*
  leave-one-out estimates. Both terms are now normalized. Since the gamma
  kernels' normalizing constant grows with the bandwidth, the old criterion
  pulled strongly towards oversmoothing. Symmetric kernels with a uniform start
  use `stats::bw.ucv()` and are unaffected.

* `bw = "ucv"` for the asymmetric kernels with a non-uniform parametric start
  took the kernel's shape from the wrong data point, so the criterion did not
  match the estimator.

* `bw = "ucv"` for the gamma, beta and Gaussian-copula kernels now searches
  around a closed-form pilot on the kernel's own bandwidth scale, avoiding
  intervals that force oversmoothing. It warns when the selected bandwidth is
  within 1% of either end of the search interval.

* `logLik()`, `AIC()` and `BIC()` for fits with a parametric start summed the
  start's densities instead of their logarithms and reported `nobs = 1`. They
  now return the start's log-likelihood and the correct AIC and BIC.

* The Gaussian-copula kernel constant is fixed, so unnormalized estimates with
  a uniform start integrate to one; normalized estimates are unchanged. Data at
  0 or 1 now give an error suggesting the beta kernel or moving the
  observations off the boundary.

* Gaussian-copula estimates at exactly 0 or 1 return their limit instead of
  `NaN`, so plots and grids that start at the boundary keep their end points.

* Fits give a clear error when the fitted start has zero, negative or
  non-finite density at an observation, instead of silently returning `NaN`.
  Numeric `bw = Inf` still returns the parametric start without this check.

* The vectorized bias-corrected beta kernel no longer reorders the values from
  its lower, interior and upper sections.

### Performance

* Built-in kernels are evaluated exactly in compiled code, so installing from
  source now needs a C++17 compiler. Gaussian and Gaussian-copula evaluation
  is much faster, and the beta and gamma kernels reuse cached data
  transformations. User-supplied kernels keep the R path.

* Compact-support kernels switch to exact sorted range sums once enough points
  have been evaluated, and the Laplace kernel uses an exact two-sweep
  recurrence, which speeds up large and repeated evaluations.

* Exact normalizing constants replace numerical integration for ordinary
  symmetric KDEs, Gaussian-copula KDEs, and Gaussian kernels with a normal
  start. This also avoids spurious normalization failures for compact
  symmetric kernels.

* `bw = "ucv"` gains compiled paths: the compiled evaluator for built-in
  kernels with a uniform start, and a log-stable pairwise path with batched
  leave-one-out estimates for normal, exponential, lognormal, gamma, Weibull
  and beta starts. Gamma and beta kernels reuse shape and normalization
  constants across observations. Invalid compiled fits fall back to the R path.

## kdensity 1.2.0 (2026-05-04)

* New closed-form bandwidth selector `"HS"`, based on Hallberg
  Szabadváry's beta reference rule (arXiv:2601.19553). It is now the
  default bandwidth for the bias-corrected `beta` kernel with a uniform
  or constant start, with a fallback heuristic when the rule's moment
  conditions are not satisfied.
* `RHE` no longer depends on the `EQL` package; the Hermite polynomials
  are evaluated inline (matches the previous numeric output to ~1e-12).
* Trimmed dependencies: `EQL` removed from `Imports`; `assertthat`,
  `testthat`, `extraDistr`, and other build-time helpers no longer
  required at install time.
* Internal cleanup: the `R/` source layout is consolidated into
  `bandwidths.R`, `kernels.R`, `starts.R`, and `helpers.R`; no
  user-visible API changes.

## kdensity 1.1.1 (2025-02-18)

* Added support for update of univariateML to version 1.5.
* There are breaking changes in the new version univariateML, but only for this package.
* The new package update is designed to work with univariateML less than 1.5 as well.

## kdensity 1.1.0 (2019-01-10)

* Journal of Open Source Software Release.
* Safer kdensity object structure.
* Cleaned up readme and documentation.

## kdensity 1.0.1 (2019-07-11)

* Fixed the `laplace` kernel.
* `kdensity`: Added error message when the normalization constant fails to be accurate enough and a tolerance parameter to control it.
* Better error messages when data is not in the support of the parametric starts.
