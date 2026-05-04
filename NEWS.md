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
