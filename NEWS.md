## kdensity 1.3.0 (2026-07-20)

### Breaking changes

* `update()` on a `kdensity` object now returns a new object and follows R's
  copy-on-modify semantics instead of modifying the caller's variable in place.
  Replace `update(kde, ...)` with `kde <- update(kde, ...)`. The `$<-` and
  `[[<-` replacement forms are unchanged.

### New features

* Added `kdensity_grid()`, an explicitly approximate linear-binning/FFT backend
  for the built-in translation-invariant kernels and the Gaussian-copula kernel
  (on the probit scale). Callable `kdensity` objects stay exact. Plotting now
  defaults to `grid = "auto"`, using verified FFT acceleration for Gaussian and
  Gaussian-copula kernels while retaining faster exact evaluation for compact
  and asymmetric kernels. Gaussian-copula grids evaluate numerically unstable
  extreme probit tails exactly. Exact grids ignore FFT-only arguments, stay
  linearly spaced, and require whole-number grid sizes.
* Added the `"GR"` gamma-reference bandwidth selector, obtained by substituting
  a fitted gamma density into Chen's AMISE-optimal bandwidth. It is not yet a
  default: the effect of final normalization and non-uniform starts still needs
  validation.

### Performance

* Built-in kernels now use an exact compiled evaluation path. Gaussian and
  Gaussian-copula evaluation is substantially faster, and beta and gamma kernels
  reuse cached data transformations. User-supplied kernels keep the R path.
* Compact-support kernels switch lazily to exact sorted-data range evaluation
  after enough cumulative work, and Laplace kernels use an exact two-sweep
  recurrence — speeding up large and repeated callable evaluations without
  adding sorting overhead to small calls.
* Exact normalization constants replace numerical integration for ordinary
  symmetric KDEs, Gaussian-copula KDEs, and Gaussian kernels with a normal
  parametric start. This also avoids spurious normalization failures for
  compact symmetric kernels.
* Unbiased cross-validation (`bw = "ucv"`) gains compiled paths: it reuses the
  compiled evaluator for built-in kernels with a uniform start, and adds a
  log-stable pairwise path for normal, exponential, lognormal, gamma, Weibull,
  and beta starts. Their leave-one-out MLEs are embedded and batched, so this
  path does not rely on unpublished `univariateML` C++ headers; non-finite or
  invalid compiled fits fall back to the existing R path.

### Bug fixes

* `bw = "ucv"` for asymmetric kernels (`gamma`, `beta`, and their bias-corrected
  variants) with a non-uniform parametric start set the kernel's shape from the
  wrong data point, making the objective inconsistent with the estimator.
  Selected bandwidths for those combinations change slightly.
* `logLik()`, `AIC()`, and `BIC()` on objects with a parametric start were
  wrong: `logLik()` summed the start's densities instead of their logarithms and
  reported `nobs = 1`. They now return the start's true log-likelihood (evaluated
  on the log scale when supported) and the correct AIC/BIC.
* Fixed the bias-corrected beta kernel's vectorized path, which previously
  reordered values from the lower, interior, and upper sections.

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
