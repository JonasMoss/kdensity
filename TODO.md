# Performance roadmap

This file is the local source of truth for performance work. Keep broad
research questions and benchmark results here. Create GitHub issues when an
item has a settled mathematical contract, a bounded implementation, and
testable acceptance criteria.

## Hardening & Result<> refactor (2026-07-20)

- [x] The compiled leave-one-out fits (`designated_loo_parameters_cpp`) no
  longer throw for data-fit failures. Each Newton fit is extracted into a
  helper returning `std::optional<...>` ("Result<>" without an error payload);
  non-convergence or degenerate data yields `std::nullopt` -> a NaN row, which
  the R layer detects (`anyNA`) and falls back on. `n < 3` returns an all-NaN
  matrix. Only the structural "unknown start id" assertion remains an
  exception. The beta Newton loop gained its missing finiteness guards here.
- [x] Enabled C++17 (`src/Makevars` / `src/Makevars.win`, `CXX_STD = CXX17`);
  R reports it as an INFO, not a NOTE.
- [x] Fixed `logLik()` / `AIC()` / `BIC()` (found while hardening tests):
  `logLik` summed the parametric start's densities instead of their logarithms
  and reported `nobs = 1`, so all three were wrong. Now returns the true
  log-likelihood. Public output values change; README regenerated. See NEWS.
- [x] Test suite hardened: new `test_degenerate.R` (graceful NaN paths) and
  `test_invariants.R` (non-negativity, closed-form normalization, reflection
  symmetry, compiled-vs-R agreement, eval-first UCV gold pin, accessor
  correctness); `test_sorted_kernels.R` gained boundary-tie coverage across
  every compact kernel.

## Open bugs (audit 2026-07-20)

Found in a review of the compiled/grid/selector dev cycle. Severity first.

- [x] **HIGH (FIXED 2026-07-20) — uniform/rectangular sorted path invents
  density at the boundary.** Window bounds now use the divided predicate;
  regression test added in `test_sorted_kernels.R`.
  `src/kernels.cpp:599-608` selected the moment-tree window with
  subtraction bounds `[lower_bound(y-h), upper_bound(y+h))`, but the direct
  path and `stats::density` include a point iff `abs((x-y)/h) <= 1`. These
  differ by one ULP at the boundary; polynomial/trig kernels vanish there so
  they are unaffected, but the uniform kernel (id 3) jumps to 0.5, so excluded
  boundary ties contribute in full. Reproducer: `x={1x20,2x60,3.5x20}`,
  `h=0.7`, `kernel="uniform"`, `y=2.7` gives 0.4286 (sorted) vs 0 (exact).
  Auto-activates at scale on tied/integer data; `test_sorted_kernels.R` misses
  it. Fix: compute the window bounds with the divided predicate.
- [x] **MEDIUM (FIXED 2026-07-20) — asymmetric-kernel non-uniform UCV
  cross-term uses swapped kernel arguments.** Now eval-first in the R fallback,
  the compiled path, and the test reference; verified against the estimator to
  ~1e-16. `R/bandwidths.R:262` and the compiled `ucv_loo_term_cpp`
  evaluate `kernel(x_j, x_i)` (shape from the other data point) while `term1`
  and the estimator use `kernel(x_i, x_j)` (shape from the eval point, Chen's
  convention). Since gamma/beta kernels are asymmetric, the LOO term is
  mis-specified, so the UCV objective is slightly wrong for the documented
  gamma+gamma / beta+beta defaults. The rewritten uniform branch is already
  eval-first, so the selector is internally inconsistent. Fix the R fallback
  and the C++ together.
- [x] LOW (FIXED 2026-07-20) — compiled leave-one-out fit aborting instead of
  degrading. `make_compiled_loo_parameters` is now wrapped in `tryCatch` and its
  result checked with `anyNA` in `bandwidths.R`; a compiled failure/NaN falls
  back to the R (univariateML) LOO path. Note: truly un-fittable data (e.g.
  `c(1e-300, 1, 2, 3, 4)` with a gamma start) still errors, but at the
  pre-existing full-data start fit that runs before bandwidth selection — that
  is out of scope and correct (the MLE does not exist). The n>=3 vs n>=2 floor
  is left as-is (n=2 non-uniform ucv is degenerate in both paths).
- [x] LOW (FIXED 2026-07-20) — beta LOO branch silent NaN. Superseded by the
  Result<> refactor above: `fit_beta_shapes_loo` now guards the non-finite
  initial value and the post-loop NaN and returns `std::nullopt`, so the R
  layer falls back via `anyNA`.
- [x] LOW (FIXED 2026-07-20) — `R/grid.R` now routes its `stats::density()`
  calls through a version-guarded `density_fft()` helper that drops
  `subdensity`/`ext` and silences the sum-of-weights warning on R < 4.4.0,
  instead of imposing a package-wide `Depends: R (>= 4.4.0)`.
- [x] LOW (FIXED 2026-07-20) — `known_normalization` now requires
  `all(x > 0 & x < 1)` for the gcopula closed form, so boundary data fails early
  (verified: constructs error at normalization instead of a broken object).
- [x] TEST (FIXED 2026-07-20) — `test_grid.R` boundary test now pins finite
  interior values against exact evaluation; the custom-kernel fallback test in
  `test_compiled_kernels.R` now checks numeric agreement, not just length/count.

## Current baseline

The first compiled pass is complete:

- All built-in kernels have an exact compiled evaluation path.
- Data-only transformations used by beta, gamma, and Gaussian-copula kernels
  are cached.
- Exact normalizers are used for ordinary symmetric KDEs, Gaussian-copula
  KDEs, and Gaussian kernels with a normal start.
- Uniform-start UCV uses compiled kernel sums.
- Normal, exponential, lognormal, gamma, Weibull, and beta starts have embedded
  batched leave-one-out MLEs and compiled log-density evaluation.
- Non-uniform UCV evaluates kernel/start ratios in one compiled, log-stable
  pairwise loop.
- Custom kernels and starts retain their R fallback.

Representative local results:

- Evaluating 5,000 observations at 1,000 points is approximately 6 times
  faster for Gaussian and Gaussian-copula kernels, 13 times faster for beta
  kernels, and 17 times faster for gamma kernels.
- Constructing the corresponding objects is approximately 9--14 times faster.
- Default gamma UCV at `n = 1,000` is approximately 5 times faster.
- In an optimized local gamma-kernel/gamma-start run at `n = 1,000`, the
  compiled non-uniform UCV term was approximately 5 times faster and the full
  bandwidth search approximately 8 times faster than the R fallback. Add
  repeated measurements to the benchmark harness before treating these ratios
  as stable.

## Completed: Parametric starts in asymmetric-kernel UCV

Parametric-start fitting and ordinary evaluation are not currently the main
bottlenecks. At `n = 1,000`, common start fits take roughly 0--1 ms, Gaussian
kernel object construction takes 1--8 ms, and evaluation at 1,000 points takes
roughly 12--13 ms across uniform, normal, Gumbel, gamma, Weibull, and lognormal
starts.

The expensive case was an asymmetric kernel with UCV. In the original
gamma-kernel/gamma-start profile at `n = 1,000`, leave-one-out fitting accounted
for only about 4% of runtime; repeated kernel and gamma-density evaluation
dominated.

The implemented contract keeps the existing `density`, `estimator`, and
`support` interface. Built-in starts may additionally carry `compiled_id`;
custom starts without it continue to use the R path.

The compiled path:

- passes a matrix of leave-one-out parameters to C++;
- evaluates the asymmetric kernel and
  `log g_{-i}(x_i) - log g_{-i}(x_j)` in the same pairwise loop.
- accumulates density ratios with log-sum-exp;
- uses sufficient statistics for normal, exponential, lognormal, gamma, and
  beta fits, and an embedded compiled Weibull fit;
- agrees with the public start estimators and corrected R criterion to test
  tolerance.

### Research fixed-start cross-validation separately

Holding the full-sample parametric start fixed during cross-validation would be
much cheaper, but changes the criterion. Do not change the default semantics
without deriving the criterion, comparing it in simulations, and exposing the
choice explicitly if it is useful.

### Analytic route for parametric starts

For a fixed start density `g`, write the target as `f = g r`. Away from the
boundary, the bias-corrected gamma kernel gives the unnormalized leading bias

`(b / 2) x g(x) r''(x)`.

After `kdensity`'s final normalization, the relevant roughness is expected to
be

`R_g = integral [x g r'' - f integral{x g r''}]^2`,

while Chen's leading integrated variance remains

`[2 sqrt(pi) n sqrt(b)]^(-1) integral x^(-1/2) f(x)`.

This suggests the non-CV bandwidth

`b = [V / {2 sqrt(pi) n R_g}]^(2/5)`.

The formula still needs a careful boundary and estimated-parameter derivation.
It also explains why simply using the fitted start as the reference is not
enough: under the exact parametric model `r = 1`, `R_g = 0`, and the
first-order optimum is infinite bandwidth.

A promising analogue of `RHE` is to expand `r` around the fitted start and
estimate a few departure coefficients:

- generalized Laguerre polynomials for a gamma start;
- Jacobi polynomials for a beta start;
- transformed Hermite polynomials for lognormal starts.

Derive `R_g` for these bases and compare the resulting selector with UCV. This
is more principled than fixed-start CV and could make non-CV selection adapt to
departures from the parametric model.

The first gamma-start derivation and simulation harness live in
`tools/gamma-start-amise.md`, `tools/gamma-start-amise.R`, and
`tools/gamma-mise-experiment.R`. A 20-repetition pilot at `n = 150` found the
Laguerre selector promising under gamma, lognormal, and Weibull targets, but
inferior to gamma-start UCV for a bimodal gamma mixture.

The bounded-support analogue now lives in `tools/beta-start-amise.md`,
`tools/beta-start-amise.R`, and `tools/beta-mise-experiment.R`. It uses a
shifted-Jacobi expansion under a fitted beta start and includes the leading
effect of final normalization. This exposed a beta-specific point: even for a
uniform start, normalization changes the leading roughness because
`integral x(1-x)f'' = f(0) + f(1) - 2`. A normalized HS comparison did not
improve a 20-repetition pilot, so there is no present reason to replace the
package's HS rule.

The beta pilot repeated the main gamma result. The analytic selector took
roughly 1--3 ms rather than 50--140 ms for UCV and was competitive for smooth,
near-beta targets, but a fixed degree-four expansion was not robust to a beta
mixture or a boundary-singular target.

The second-stage work is summarized in
`tools/series-selector-findings.md`. Completed research infrastructure now
includes:

- fixed oracle-MISE grids for uniform and parametric starts;
- degree 2, 4, 6, and 8 comparisons;
- hard, soft, nonnegative-garrote, global-test, and AIC/BIC rules;
- an analytic composite-null covariance correction that projects coefficient
  noise off the fitted gamma or beta score space;
- an explicit diagnosis of the beta boundary square-integrability problem;
- experimental shrinking-level global gates and a beta boundary guard that
  switches to uniform-start HS/UCV;
- parametric-bootstrap calibration of the composite-null global statistic,
  with exact-null comparisons at `n = 75, 150, 300` and added close/interior
  mixture scenarios;
- a continuous bootstrap-tail transition
  `w = (T / c)^p / {1 + (T / c)^p}`, where `c` is the fitted-null bootstrap
  critical value;
- gamma inverse-roughness interpolation,
  `b_smooth = b_Laguerre / w^(2/5)`, and beta density averaging between the
  parametric fit and finite-bandwidth Jacobi estimate;
- small smooth-transition pilots at `n = 75, 150, 300`, using a gentle
  `p = 1` gamma transition and a sharper `p = 4` beta transition;
- transformed midpoint MISE quadrature whose boundary results are stable as
  the evaluation grid is refined.

Higher degrees close much of the oracle gap on mixtures, but no automatic
selection rule dominates smooth and mixture cases. The beta boundary guard is
a strong small-pilot improvement, while the shrinking-level gate deliberately
trades mixture power for protection of exact models. The continuous
transition recovers much of that lost power without the rare catastrophic
null expansions of a hard gate, but still trails always-on degree-eight
smoothing on strong alternatives. Keep these selectors out of the package API
for now. Remaining work, if this line is pursued further:

- try blockwise risk estimation rather than coefficientwise shrinkage;
- include more mixtures, monotone densities, and heavy tails;
- finish a uniform boundary expansion for the parametric-start estimators;
- calibrate the beta shape and endpoint-power guardrails on fresh scenarios;
- replace per-fit bootstrap with tabulated critical values or a response
  surface indexed by sample size, degree, and fitted start parameters;
- extend the smooth-transition comparison:
  - bag the binary decision over bootstrap samples and use the selection
    frequency as a mixing weight;
  - replace the beta shape and endpoint cutoffs by smooth transition weights,
    while retaining a hard safety stop when the fitted start is not square
    integrable;
  - compare smooth weights based on calibrated p-values, excess over a
    bootstrap critical value, and estimated risk;
- run several sample sizes with at least 1,000 repetitions before
  reconsidering package exposure.

## P1: More exact kernel algorithms

Completed in `src/kernels.cpp`, with benchmarks in
`tools/exact-kernel-benchmark.md`:

- rectangular, triangular, Epanechnikov, biweight, triweight, tricube, cosine,
  and optimal-cosine kernels use exact sorted range queries over cached local
  moments;
- the Laplace kernel uses exact forward and backward exponential recurrences;
- evaluation switches lazily after cumulative work reaches an empirical
  break-even threshold, so small and scalar calls retain direct summation;
- the local-moment tree remains stable when the data have a large location
  offset, unlike global raw-prefix moments.

## P1: Fast grid evaluation

Completed in `R/grid.R`, with benchmarks in `tools/grid-benchmark.md`:

- `kdensity_grid()` provides an explicitly approximate linear-binning/FFT
  backend for the seven translation-invariant kernels supported by
  `stats::density()`;
- `plot`, `lines`, and `points` opt in through `grid = "fft"` or use
  `grid = "auto"` for an exact fallback; callable `kdensity` objects remain
  exact;
- the Gaussian-copula kernel uses its exact Gaussian-convolution
  representation on the probit scale;
- `fft_n` controls binning resolution, `truncate` controls Gaussian padding,
  and accuracy tests compare every supported path with exact evaluation.

Beta and gamma kernels are not translation-invariant because their parameters
change with the evaluation point. They do not have a direct ordinary-FFT
convolution path. Keep exact compiled summation as the default while
investigating controlled binning or exponential-sum approximations separately.

## P1: Bandwidth selection

- The explicit `GR` selector now substitutes a fitted gamma reference into
  Chen's AMISE-optimal bandwidth for the bias-corrected gamma kernel. Validate
  how the leading bias changes under `kdensity`'s final normalization and under
  non-uniform parametric starts before considering it as a default.
- Compare `GR` with UCV across monotone, heavy-tailed, multimodal, and
  near-boundary densities.
- Preserve UCV as an explicit robust option.
- Keep the existing closed-form beta selector and expand its benchmark
  coverage.

## Benchmark and acceptance harness

- Measure start fitting, bandwidth selection, normalization, construction, and
  evaluation separately.
- Cover `n = 100, 1,000, 10,000` and both scalar and 1,000-point evaluation.
- Benchmark every built-in kernel with uniform and representative compatible
  parametric starts.
- Record maximum absolute and relative error against the R reference.
- Treat speedups as acceptable only when tests, examples, and `R CMD check`
  remain clean.

## Candidate GitHub issues

Create these when work begins, linking back to this roadmap:

1. Add a reproducible performance and accuracy harness.
2. ~~Extend start objects with optional compiled/log-density capabilities.~~
3. ~~Compile non-uniform asymmetric-kernel UCV for designated starts.~~
4. ~~Add exact fast algorithms for compact and Laplace kernels.~~
5. ~~Add an optional FFT/binning grid backend.~~
6. Research and validate a gamma-reference bandwidth selector.
