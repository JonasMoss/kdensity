# Performance roadmap

This file tracks unfinished performance and selector work. Completed changes
belong in `NEWS.md` and the git history. Create GitHub issues when an item has a
settled mathematical contract, a bounded implementation, and testable
acceptance criteria.

## Release status

There are no known blockers for the current CRAN release. The compiled kernel
paths, leave-one-out UCV implementation, exact sorted evaluators, optional FFT
grid backend, and explicit `GR` selector are part of the tested package.

The research below is not release-ready and is not exposed through the package
API unless stated otherwise.

## P1: Reproducible performance and accuracy harness

- [ ] Measure start fitting, bandwidth selection, normalization, construction,
  and evaluation separately.
- [ ] Cover `n = 100, 1,000, 10,000` and both scalar and 1,000-point
  evaluation.
- [ ] Benchmark every built-in kernel with uniform and representative
  compatible parametric starts.
- [ ] Use repeated measurements and report uncertainty rather than isolated
  speedup ratios.
- [ ] Record maximum absolute and relative error against the exact or R
  reference.
- [ ] Keep benchmark acceptance conditional on 100% coverage, clean tests,
  examples, and `R CMD check`.

Existing starting points:

- `tools/verification/fast-grids/verify.R`
- `tools/exact-kernel-benchmark.R`
- `tools/grid-benchmark.R`
- `tools/gamma-mise-experiment.R`
- `tools/beta-mise-experiment.R`

## P1: Validate the gamma-reference bandwidth selector

`GR` is release-ready as an explicit opt-in selector, but not as a replacement
for an existing default.

- [ ] Finish the leading-bias derivation under final normalization and
  non-uniform parametric starts.
- [ ] Compare `GR` with UCV across monotone, heavy-tailed, multimodal, and
  near-boundary densities.
- [ ] Run repeated simulations across several sample sizes.
- [ ] Keep UCV available as the robust explicit option.
- [ ] Consider changing a default only if the derivation and simulations show
  a clear, stable advantage.

## Research only: Analytic parametric-start selectors

The experimental Laguerre and shifted-Jacobi selectors are substantially
faster than UCV and can work well for smooth near-parametric targets, but no
automatic rule currently performs reliably across exact models, mixtures, and
boundary-singular targets.

Before considering package exposure:

- [ ] Try blockwise risk estimation rather than coefficientwise shrinkage.
- [ ] Include more mixtures, monotone densities, and heavy tails.
- [ ] Finish a uniform boundary expansion for the parametric-start estimators.
- [ ] Calibrate beta shape and endpoint-power guardrails on fresh scenarios.
- [ ] Replace per-fit bootstrap calibration with validated tabulated critical
  values or a response surface.
- [ ] Compare smooth transition weights based on bootstrap selection
  frequencies, calibrated p-values, critical-value excess, and estimated risk.
- [ ] Retain a hard safety stop when a fitted beta start is not square
  integrable.
- [ ] Run several sample sizes with at least 1,000 repetitions before
  reconsidering package exposure.

Research notes and implementations:

- `tools/gamma-start-amise.md`
- `tools/beta-start-amise.md`
- `tools/series-selector-findings.md`

## Research only: Faster asymmetric-kernel grids

Beta and gamma kernels are not translation-invariant because their parameters
change with the evaluation point, so they do not have an ordinary FFT
convolution path.

- [ ] Investigate controlled binning or exponential-sum approximations.
- [x] Define per-call accuracy, fallback, support-boundary, and release
  contracts in `tools/verification/asymmetric-grids/CONTRACT.md`.
- [ ] Prototype the one-parameter gamma and beta exponential sums behind the
  contract; keep corrected-beta boundary regions exact initially.
- [ ] Build the asymmetric-grid accuracy and performance verifier described by
  the contract before adding a package backend.
- [ ] Keep exact compiled summation as the default.

## Candidate GitHub issues

Create these when work begins:

1. Add a reproducible performance and accuracy harness.
2. Validate the gamma-reference bandwidth selector across target families.
