# Approximate asymmetric-grid contracts

Status: research contract. The exact compiled evaluator remains the default.

This document separates two questions that should not be conflated:

1. What must one approximate call guarantee?
2. What evidence is required before an approximate backend is exposed or used
   automatically?

An empirical comparison against a few exact grids is release evidence, but it
is not the per-call accuracy guarantee.

## 1. Exact semantic target

For observations `x_i`, fitted parametric-start density `r`, fitted
normalization constant `Z`, and an asymmetric kernel density `q_y,h`, the
package evaluates

```text
f(y) = r(y) / Z * sum_i q_y,h(x_i) / (n * r(x_i)).
```

The approximation may reorganize this positive weighted sum. It must not
change the kernel shapes, refit or discretely renormalize the result, alter the
parametric-start correction, or replace `Z`.

The reference branch rules are the rules in `src/kernels.cpp`:

- `gamma_biased`: gamma shape `a(y) = y / h + 1` and scale `h`;
- `gamma`: `a(y) = (y / h)^2 / 4 + 1` when `y < 2h`, otherwise
  `a(y) = y / h` (Chen's modified kernel; corrected in 2.0.0, earlier
  versions used `a(y) + 1`);
- `beta_biased`: beta shapes `a(y) = y / h + 1` and
  `b(y) = (1 - y) / h + 1`;
- `beta`: use the lower corrected shapes when `y < 2h`, the interior
  shapes when `2h <= y < 1 - 2h`, and the upper corrected shapes otherwise.

That branch precedence matters when the two beta cut points meet or cross.
The compiled evaluator, rather than a mathematically similar reformulation,
is the release-test oracle.

## 2. Per-call approximation contract

Let `f_hat[j]` be an approximate value and `f[j]` the exact semantic target at
the same returned grid point. For requested non-negative tolerances `atol` and
`rtol`, at least one of which is positive, every approximated point must obey

```text
abs(f_hat[j] - f[j]) <= atol + rtol * abs(f[j]).
```

This is a mixed bound because a pure relative bound is undefined at exact
zeros and unhelpful in underflow-level tails.

### Required certificate

The backend must construct, without evaluating the full exact sum:

- `E[j]`, an upper bound on the absolute approximation error; and
- `L[j]`, a non-negative lower bound on the exact density.

An approximate point is accepted only when

```text
E[j] <= atol + rtol * L[j].
```

Since `L[j] <= f[j]`, this condition proves the public mixed bound. Bounds may
be shared by a query block, but they must hold for every point in that block.
If a bound is unavailable, non-finite, or too large, the point or block must
fall back to exact compiled summation.

The error budget must include every approximation layer that is used:

- source grouping or binning;
- interpolation or truncated series;
- transform truncation and evaluation-grid interpolation; and
- floating-point evaluation and summation error.

Tuning a bin count or expansion order on the verification scenarios is not a
certificate. The chosen order may be data-dependent, but its bound must be
evaluated for the current fit, bandwidth, and query block.

### Result invariants

An approximate or hybrid call must also:

- return the same grid, length, order, and names as exact evaluation;
- return finite, non-negative values wherever the exact result is finite and
  non-negative;
- preserve exact zeros implied by support;
- be deterministic for fixed inputs and documented thread settings; and
- respond to user interrupts with the same granularity expected of the exact
  compiled path.

Negative approximations must not be silently hidden. A value may be clamped
to zero only if the clamp displacement is added to `E[j]` and the mixed bound
still passes; otherwise it falls back to exact evaluation.

### Boundary and exceptional inputs

The first implementation should use exact summation for:

- beta observations at `0` or `1` and gamma observations at `0`, unless their
  contribution has a separately derived endpoint formula;
- beta query endpoints `0` and `1`, gamma query point `0`, and the exact branch
  knots `2h` and `1 - 2h`;
- corrected-beta boundary regions until their two-parameter approximation has
  its own certificate;
- non-finite transformed coordinates or unsafe exponent ranges; and
- points where the fitted start multiplier or normalization cannot safely
  translate a correction-sum bound into a density bound.

The automatic route must use exact evaluation outside the kernel support or
when the exact evaluator itself has non-finite semantics. Approximation must
not turn an exact `Inf`, `NaN`, or domain error into an ordinary finite value.

Exact boundary evaluation is an initial implementation rule, not a claim that
the boundary formulas are intrinsically unapproximable. A boundary branch can
be promoted independently after it satisfies this same contract.

## 3. Candidate computational reductions

The package formulas expose useful exponential sums, despite not being
ordinary convolutions.

For either gamma kernel, with `t = a(y) - 1`, the varying part is

```text
sum_i w_i * exp(t * log(x_i) - x_i / h),
```

times a known function of `t` and `h`. Both the linear and corrected gamma
shape paths therefore query the same one-parameter exponential sum.

For `beta_biased`, put `t = y / h`. After taking out known beta-function and
bandwidth factors, the varying part is

```text
sum_i w_i * (1 - x_i)^(1 / h) * exp(t * logit(x_i)).
```

The interior branch of corrected `beta` has the same one-parameter structure
with different fixed powers. Its corrected lower and upper branches trace a
two-parameter curve in `(a - 1, b - 1)` and need either a two-dimensional
low-rank scheme or exact fallback.

The first prototypes worth comparing are:

1. an unequally spaced fast Laplace/exponential-sum evaluator;
2. an adaptive source/query tree with bounded Taylor or Chebyshev blocks; and
3. controlled transformed-coordinate binning with an analytic bin remainder.

Ordinary linear binning plus FFT is not a candidate on the original scale:
these kernels change shape with the evaluation point. Any prototype that only
reports an observed error, without a current-call remainder bound, stays in
the research harness.

## 4. Proposed tolerance profile

The core contract is parameterized by explicit `atol` and `rtol`. A first
automatic profile should be studied at

```text
rtol = 1e-4
atol = 1e-10 * scale
```

where `scale` is the maximum finite density from a small, fixed set of exactly
evaluated pilot grid points. If no positive finite scale is available, the
whole call is exact. The pilot set and its cost must be included in timings.

The research harness must also exercise `rtol` values `1e-3`, `1e-5`, and
`1e-7`. The default profile is not frozen until the speed/accuracy curves show
that it is useful. The initial package API need not expose tolerance controls;
internal controls are sufficient to verify the contract without adding public
surface area.

## 5. Release-validation contract

The verifier must be reproducible, live under this build-ignored directory,
write machine-readable case results, and exit nonzero on a promotion failure.
It must record the requested tolerance, certified bound, observed error,
fraction of exact points, preparation time, evaluation time, and peak memory
where available.

### Accuracy matrix

Cover all four built-in asymmetric kernels across:

- sample sizes `32`, `100`, `1,000`, and `10,000`;
- output grids of `64`, `512`, and `2,048` points;
- narrow, selected, and wide bandwidths, including beta bandwidths for which
  the two correction regions nearly meet or cross;
- full-support and interior subranges, with grids that contain endpoints and
  branch knots exactly;
- representative built-in starts and at least one custom positive start;
- gamma targets concentrated near zero, exponential, lognormal/heavy-tailed,
  multimodal, tied, and containing exact zero observations;
- beta targets that are interior, skewed, U-shaped, boundary-heavy, tied, and
  contain exact `0` or `1` observations; and
- start-correction weights spanning ordinary and severe dynamic ranges.

Use fixed adversarial cases plus seeded randomized cases. Compare against the
compiled exact evaluator in every case and against an independent scalar
`dbeta`/`dgamma` implementation on a smaller golden subset. A high-precision
reference may be used by the developer harness, but must not become a package
dependency.

Promotion requires:

- zero violations of the requested mixed pointwise tolerance;
- zero cases where observed error exceeds the recorded certificate, allowing
  only an explicitly derived reference-roundoff allowance;
- exact agreement on routing and support/exception behavior;
- finite and non-negative invariant checks to pass in every applicable case;
- relative integrated absolute error and relative mass error no worse than
  the requested pointwise profile permits; and
- package test coverage to remain 100% for all reachable package code.

### Performance matrix

Measure warmed exact and hybrid calls separately from cache construction, and
also measure total first-call cost. Use at least seven repetitions and report
medians and a robust spread. The pilot, certificate construction, and exact
fallback work count toward the hybrid time.

Automatic promotion is per kernel and requires, for `n >= 1,000` and a
512-point grid:

- at least a two-fold geometric-mean total speedup across representative
  scenarios;
- a 10th-percentile speedup of at least `1.25`;
- no adversarial case selected by `"auto"` to be slower than exact by more
  than 10%; and
- small calls routed to exact unless the total approximate call remains below
  10 ms and wins reliably.

The routing decision must use information available before doing most of the
approximate work. A prototype that is fast only after excluding its exact
fallback or preprocessing cost is not ready for `"auto"`.

## 6. Promotion stages

1. **Research only:** scripts under `tools/verification/asymmetric-grids/`;
   no package path and no CRAN-visible change.
2. **Explicit hybrid:** an opt-in approximate method may be exposed after the
   per-call and release-accuracy contracts pass on supported branches. All
   unsupported branches remain exact.
3. **Automatic per kernel:** `method = "auto"` may select a hybrid only after
   that kernel also passes the performance contract on supported platforms.
4. **Broader branch promotion:** corrected-beta boundary branches are promoted
   separately; they do not block useful one-parameter paths.

Any numerical bug that violates the pointwise contract removes the affected
kernel or branch from automatic routing until a new verification report
passes. Exact compiled summation remains the permanent fallback and semantic
oracle.

## Primary references

- Chen, S. X. (1999), [Beta kernel estimators for density
  functions](https://songxichen.com/Uploads/Files/Publication/Chen-CSD-99.pdf).
- Chen, S. X. (2000), [Probability density function estimation using gamma
  kernels](https://doi.org/10.1023/A:1004165218295).
- Lee, D. and Gray, A. G. (2012), [Faster Gaussian Summation: Theory and
  Experiment](https://arxiv.org/abs/1206.6857). The reusable idea here is the
  combination of positive weights, query lower bounds, and accumulated
  absolute-error budgets, not its Gaussian-specific expansion.
- Andersson, F. (2013), [Algorithms for unequally spaced fast Laplace
  transforms](https://doi.org/10.1016/j.acha.2012.11.005).
- Blanchard, P., Higham, N. J., and Mary, T. (2020), [A Class of Fast and
  Accurate Summation Algorithms](https://doi.org/10.1137/19M1257780).
