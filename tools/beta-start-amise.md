# Beta-start AMISE research note

This note mirrors the gamma-start experiment for the bias-corrected beta
kernel. It is research code, not a package selector. The implementation is in
`beta-start-amise.R`; the paired MISE harness is in
`beta-mise-experiment.R`.

The starting point is Chen's second beta-kernel estimator. Away from the two
boundary strips, its kernel at the evaluation point `x` is a beta density with
parameters `x / b` and `(1 - x) / b`. Its mean is exactly `x`, and its
variance is

```text
b x(1-x) / (1+b) = b x(1-x) + O(b^2).
```

Consequently, Chen's leading interior bias is

```text
(b / 2) x(1-x) f''(x),
```

and its leading integrated variance is

```text
V / (2 sqrt(pi) n sqrt(b)),

V = integral [x(1-x)]^(-1/2) f(x) dx.
```

Chen also shows that the special boundary strips do not change those leading
integrated terms. See [Chen
(1999)](https://doi.org/10.1016/S0167-9473(99)00010-9). The current `HS`
selector substitutes a beta reference density into this AMISE; see [Hallberg
Szabadvary
(2026)](https://doi.org/10.1080/10618600.2026.2670662).

## 1. A beta parametric start

Let `g` be a beta start and write the target as

```text
f(x) = g(x) r(x).
```

The raw parametrically guided estimator smooths the density ratio and
multiplies it by `g`. Its leading interior expectation is therefore

```text
g(x) {r(x) + (b / 2) x(1-x) r''(x)}.
```

Define

```text
D(x) = x(1-x) g(x) r''(x),
A    = integral D(x) dx.
```

Before final normalization the leading bias is `(b / 2) D`. The raw estimate
has total mass `1 + (b / 2) A + o(b)`, so `kdensity`'s final normalization
changes the leading bias to

```text
(b / 2) {D(x) - A f(x)}.
```

The relevant roughness is thus

```text
J_g = integral {D(x) - A f(x)}^2 dx.
```

The local start ratios cancel in the leading variance term, leaving Chen's
`V`. Minimizing

```text
AMISE(b) = (b^2 / 4) J_g + V / (2 sqrt(pi) n sqrt(b))
```

gives

```text
b* = {V / (2 sqrt(pi) n J_g)}^(2/5).
```

This is the same structural result as for a gamma start, with `x` replaced by
`x(1-x)`.

## 2. Why normalization matters more than for gamma

For a uniform start, `D = x(1-x) f''`. Integration by parts gives

```text
A = f(0) + f(1) - 2.
```

For a smooth beta reference with both shapes greater than `3/2`, the endpoint
densities vanish and `A = -2`. Thus the package's normalized estimator has
leading bias proportional to

```text
x(1-x) f''(x) + 2 f(x),
```

whereas the published HS rule uses `x(1-x) f''(x)`. The helper
`normalized_beta_reference_bandwidth()` computes this alternative exactly,
using the same method-of-moments beta fit as `compute_hs_bandwidth()`.

This correction is mathematically real, but the small pilot below does not
show a finite-sample improvement from changing HS. It should remain an
experimental comparison.

## 3. Shifted Jacobi expansion

For

```text
g(x) = x^(a-1) (1-x)^(c-1) / B(a,c),
```

the natural orthogonal basis is

```text
phi_k(x) =
  P_k^(c-1,a-1)(2x-1) / sqrt(H_k),
```

where `P` is a Jacobi polynomial and

```text
H_k =
  Gamma(k+a) Gamma(k+c) Gamma(a+c)
  --------------------------------------------------------------- .
  (2k+a+c-1) Gamma(k+1) Gamma(k+a+c-1) Gamma(a) Gamma(c)
```

These functions are orthonormal under `g`. Expand

```text
r_K(x) = sum_(k=0)^K theta_k phi_k(x),
theta_k = E_f[phi_k(X)].
```

The empirical coefficients are just sample means. The second derivative is
cheap because

```text
d^2/dx^2 P_k^(alpha,beta)(2x-1)
  = (k+alpha+beta+1)(k+alpha+beta+2)
      P_(k-2)^(alpha+2,beta+2)(2x-1).
```

Write `R` and `S` for the ordinary polynomial coefficients of `r_K` and
`r_K''`. Then

```text
A = B(a,c)^(-1) integral x^a (1-x)^c S(x) dx
```

and

```text
D - A f =
  x^(a-1) (1-x)^(c-1) / B(a,c)
    {x(1-x) S(x) - A R(x)}.
```

The expression in braces is a finite polynomial `Q(x)`. If
`Q(x)^2 = sum q_j x^j`, the roughness is the finite sum

```text
J_g = B(a,c)^(-2) sum_j q_j B(2a-1+j, 2c-1).
```

No numerical differentiation, quadrature, or bandwidth optimization is
needed. The implementation requires fitted shapes greater than `1/2`, which
ensures the normalized squared-density integrals exist.

## 4. Degenerate and near-parametric cases

Under an exact beta model, `r = 1`, hence `r'' = 0`, `J_g = 0`, and the
first-order AMISE optimum is infinite bandwidth: use the parametric start
without smoothing. Sampling noise makes the estimated higher-order Jacobi
coefficients nonzero, so an unthresholded finite expansion can instead return
a large finite bandwidth.

The experiment includes:

- degree two and degree `K` expansions;
- a hard threshold at `1 / sqrt(n)` for coefficients of degree two and above;
- normalized and raw roughness versions;
- the exact parametric fit, HS, and uniform- and beta-start UCV.

The corrected beta boundary construction also requires `b < 1/4`, because its
lower and upper boundary strips must not overlap. The experiment caps finite
Jacobi bandwidths just below `1/4` and records every cap. An infinite
bandwidth remains the exact parametric fit.

Fitting `g` and estimating the expansion on the same data is harmless to the
fixed-alternative plug-in calculation at first order, but it matters in the
near-parametric regime. A full theory should project coefficient noise away
from fitted beta score directions or use an explicit risk estimate. That is
one reason not to expose this as a production selector yet.

## 5. Reproducible MISE experiment

Run, for example,

```sh
Rscript tools/beta-mise-experiment.R \
  --reps=1000 --n=150 --degree=4 --points=1000
```

The default scenarios are:

- `Beta(3,5)`;
- boundary-singular `Beta(0.8,0.8)`;
- a two-component beta mixture;
- a Kumaraswamy alternative;
- a logistic-normal alternative.

Each method receives the same sample within a scenario and replication. Raw
and summarized CSV files go to `tools/results/`, which is ignored by Git.
Selectors, scenarios, and output paths can be restricted with command-line
options; see the defaults in the script.

The evaluation grid is sine-squared to resolve both boundaries. All default
targets have square-integrable densities, but ISE estimates for the
boundary-singular case are intrinsically more numerically delicate than those
for the smooth cases.

## 6. Checks and preliminary results

The implementation has been checked in three ways:

1. shifted Jacobi squared norms agree with numerical integration;
2. analytic `A` and `J_g` agree with direct numerical quadrature;
3. the raw beta-reference helper agrees with the existing closed-form HS
   bandwidth to numerical tolerance.

A 20-repetition smoke pilot at `n = 150`, degree four, and 500 grid points
gave the following MISE values. These are diagnostic, not publication-quality
Monte Carlo estimates.

| Target | Leading methods by MISE |
|---|---|
| Beta(3,5) | parametric 0.01364; thresholded Jacobi-4 0.01654; beta-start UCV 0.01935; normalized Jacobi-4 0.01978; HS 0.02071 |
| Beta(0.8,0.8) | uniform UCV 0.06260; HS fallback 0.07678; beta-parametric and Jacobi methods about 0.279--0.287 |
| Beta mixture | HS 0.05035; raw Jacobi-4 0.06254; normalized/thresholded Jacobi-4 about 0.0657; uniform UCV 0.14845 |
| Kumaraswamy | parametric 0.01238; thresholded Jacobi-4 0.01477; beta-start UCV/Jacobi-2 about 0.01563; HS 0.02051 |
| Logistic-normal | parametric 0.01705; beta-start UCV/Jacobi-2 about 0.01913; thresholded Jacobi-4 0.01999; HS 0.02406 |

The normalization-adjusted HS comparison gave MISE 0.02095 under `Beta(3,5)`
and 0.02435 under Kumaraswamy, both slightly worse than ordinary HS. Its shape
condition also failed for all boundary and mixture samples and half of the
logistic-normal samples. There is currently no empirical reason to replace
HS.

Jacobi selection and fitting generally took about 1--3 ms. Beta-start UCV
took about 80--140 ms, and uniform UCV about 50--75 ms. The beta result
therefore repeats the gamma pattern:

- the analytic parametric-start route is fast and competitive for smooth,
  near-reference targets;
- fixed low-degree expansions do not adapt reliably to mixtures or boundary
  singularities;
- hard thresholding helps near the parametric model but is not uniformly
  best;
- UCV is still useful as a robust fallback, particularly at boundaries.

The shared continuation is recorded in `series-selector-findings.md`. It adds
oracle bandwidth grids, degrees 2--8, hard/soft/garrote rules, an analytic
fitted-score covariance correction, global smooth tests, and AIC/BIC degree
selection.

Degree 6--8 closes most of the beta-start oracle gap on the beta mixture, but
the corresponding noise under exact and near-beta targets remains difficult
to select away uniformly. The beta boundary case also has a separate
square-integrability problem when fitted shapes approach `1/2`. The Jacobi
selector should remain a research tool.

The experiment now makes that fallback concrete in
`beta_JacobiK_guarded`. It declines the beta start when either fitted shape is
at most 1.1, uses an order-statistic endpoint-power estimate to distinguish
evidence of a singular endpoint, and chooses uniform-start UCV for that case
or HS otherwise. In a fresh 500-replication, `n = 150` run it guarded 499/500
`Beta(0.8,0.8)` samples and all 500 beta mixtures, with MISE 0.0390 and
0.0465. It did not activate on the exact-beta, Kumaraswamy, or logistic-normal
samples. The experiment now uses transformed midpoint quadrature; the
boundary MISE is stable from 300 through 1,000 points, unlike the old endpoint
trapezoids. The exact cutoffs still need broader scenarios and sample sizes
before package exposure.

The bootstrap-calibrated guarded method reduces exact-beta false expansions:
the counts over 300 repetitions were 4, 1, and 0 at `n = 75, 150, 300`,
compared with 10, 10, and 12 for the asymptotic gate. A new interior beta
mixture shows the cost of making that calibration a hard switch. Bootstrap
selection rose from 40% to 86% to 100% over those sample sizes, and MISE was
0.375, 0.117, and 0.028; always using degree eight gave 0.090, 0.045, and
0.028. This strongly motivates smooth model averaging while keeping the hard
boundary safety guard.
