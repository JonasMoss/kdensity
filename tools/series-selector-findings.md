# Laguerre/Jacobi selector findings

This note records the shared second-stage work on the experimental
gamma-start Laguerre and beta-start Jacobi bandwidth selectors. The
distribution-specific AMISE derivations remain in `gamma-start-amise.md` and
`beta-start-amise.md`.

## 1. Oracle-MISE grids

Both experiment scripts accept `--oracle_points`. A positive value of at least
two evaluates fixed logarithmic bandwidth grids for every replication under:

- the uniform start;
- the matching gamma or beta parametric start;
- an additional infinite-bandwidth candidate for the parametric start.

The scripts write raw grid results, per-bandwidth MISE summaries, and selected
oracle-MISE optima. These are fixed bandwidths minimizing average ISE across
replications, not the more optimistic sample-wise ISE minimizers.

For example:

```sh
Rscript tools/gamma-mise-experiment.R \
  --reps=1000 --n=150 --points=1000 --oracle_points=31

Rscript tools/beta-mise-experiment.R \
  --reps=1000 --n=150 --points=1000 --oracle_points=31
```

A 20-replication pilot at `n = 150` gave:

| Target | Estimator family | Oracle bandwidth | Oracle MISE |
|---|---:|---:|---:|
| Gamma(3,1.5) | gamma start | Inf | 0.00208 |
| Gamma mixture | gamma start | 0.150 | 0.00354 |
| Lognormal | gamma start | 0.105 | 0.00382 |
| Weibull | gamma start | Inf | 0.00285 |
| Beta(3,5) | beta start | Inf | 0.01364 |
| Beta mixture | beta start | 0.0162 | 0.05576 |
| Kumaraswamy | beta start | Inf | 0.01238 |
| Logistic-normal | beta start | Inf | 0.01705 |

The finite gamma boundary optimum reached the upper edge of the pilot grid,
so that one number is intentionally not reported as settled.

These results separate two problems. The parametric-start estimator family is
useful: it has a substantially better oracle for the gamma and beta mixtures.
Most of the remaining gap is bandwidth-selector error, not an inherent failure
of parametric-start smoothing.

## 2. Estimated-start coefficient covariance

The original threshold `1 / sqrt(n)` treats the orthogonal basis as though the
start parameters were fixed. They are estimated from the same sample. Let
`s(X)` be the parametric score, `I` its Fisher information, and `phi` the
vector of nonconstant Laguerre or Jacobi functions. Under the exact
parametric model,

```text
C = E_g[phi(X) s(X)^T],

sqrt(n) theta_hat
  -> Normal(0, Omega),

Omega = I_basis - C I^(-1) C^T.
```

The subtraction is the projection of the smooth departures onto the fitted
parametric score space. This is the usual composite-null correction for
Neyman smooth tests; compare [Rayner and Best
(1986)](https://doi.org/10.1093/biomet/73.2.437).

All entries of `C` and `I` are analytic here:

- gamma uses gamma moments, digamma differences, and the gamma
  shape/rate Fisher information;
- beta uses beta moments, digamma differences, and the beta Fisher
  information.

The selectors now expose `coefficient_sd` and
`coefficient_covariance`, and `threshold_scale = "fitted"` uses their
diagonal rather than one. In 500 simulations at `n = 400`, the empirical and
predicted standard deviations of coefficients two through four were:

| Start | Empirical | Projected analytic |
|---|---|---|
| Gamma | 0.546, 0.884, 0.900 | 0.568, 0.906, 0.966 |
| Beta | 0.249, 0.968, 0.988 | 0.252, 0.987, 0.982 |

At degree eight and `n = 150`, a quadratic-form five-percent test rejected an
exact gamma start in 3.6% of 1,000 runs and an exact beta start in 6.1%.
That is sufficiently close for research comparisons, although finite-sample
calibration would be needed for a public selector.

## 3. Shrinkage and degree selection

The helpers support four coefficient rules:

- none;
- hard thresholding;
- soft thresholding;
- positive-part nonnegative-garrote shrinkage.

The experiments also include:

- `*_adaptive`: fitted-score hard thresholding with a simultaneous
  `sqrt(2 log K)` cutoff;
- `*_global`: a five-percent quadratic smooth test choosing between the
  parametric fit and the full requested degree;
- `*_bic` and `*_aic`: data-driven smooth-test order selection, including
  degree zero as the parametric fit.

The AIC/BIC order criterion follows the structure of data-driven Neyman tests:
for each prefix it maximizes a score quadratic minus either `2` or `log(n)`
per effective component. See [Ledwina
(1994)](https://doi.org/10.1080/01621459.1994.10476834).

## 4. Degree pilot

A paired 50-replication pilot compared degrees 2, 4, 6, and 8 at `n = 150`.
The output is not committed; it can be recreated with the experiment scripts
and the same seed.

For gamma:

- exact gamma favored the parametric fit (MISE 0.00192); degree-two soft
  shrinkage was close at 0.00194;
- unshrunk normalized Laguerre MISE on the gamma mixture improved from
  0.01173 at degree two to 0.00635 at degree eight, versus oracle 0.00354;
- lognormal MISE improved from roughly 0.0081 at degree two to 0.00426 at
  degree eight, versus oracle 0.00382;
- Weibull still favored the parametric fit;
- the degree-eight global rule protected the exact gamma and Weibull cases,
  but selected smoothing in only 58% of mixture samples and had MISE 0.00794.

For beta:

- exact beta favored the parametric fit (MISE 0.01317);
- unshrunk Jacobi MISE on the beta mixture fell from 0.23291 at degree two to
  0.05418 at degree eight; omitting the normalization correction gave 0.05321,
  versus beta-start oracle 0.05576 in the smaller oracle pilot;
- degree-eight global selection detected all 50 mixture samples and attained
  MISE 0.05418, while selecting smoothing in only 4--8% of the smooth
  near-parametric scenarios;
- BIC protected the parametric cases but often selected too low a degree for
  the mixture; AIC reached 0.05803 on the mixture but was noisier under the
  smooth targets.

The main conclusion is not “use degree eight.” Increasing degree reduces
truncation error on mixtures but increases roughness noise under smooth
targets. Simple soft, hard, AIC, BIC, and global rules move along that
tradeoff; none dominates in this pilot.

## 5. Two explicit guardrails

The two obvious failure modes can be handled more directly than a single
universal selector suggests.

First, the experiment methods `gamma_LaguerreK_conservative` and
`beta_JacobiK_guarded` use a global-test level

```text
alpha_n = min(0.01, 1 / n).
```

Thus the probability of smoothing an exact parametric model is designed to
vanish with sample size. In a fresh 500-replication run at `n = 150`, it
selected 13 exact-gamma, 6 boundary-gamma, 103 gamma-mixture, 208 lognormal,
and 3 Weibull samples. MISE was 0.00174, 0.01744, 0.00963, 0.00759, and
0.00334, respectively. The 2.6% exact-gamma selection rate is appreciably
larger than the nominal 0.67% level: shrinking `alpha` reduces false
expansions, but does not repair the far tail of the asymptotic chi-squared
calibration. Its lower mixture power also costs MISE relative to always using
degree eight (0.00635 in the earlier paired pilot). This is a policy knob, not
a uniformly free improvement; a production test would need parametric-
bootstrap or separately tabulated critical values.

The experiment now includes that first option as
`gamma_LaguerreK_bootstrap` and `beta_JacobiK_bootstrap_guarded`. Each method
simulates from the fitted start, refits the start and series covariance in
every bootstrap sample, and uses the finite-sample global-statistic tail. The
helper preserves the R RNG state so adding a bootstrap method does not change
subsequent Monte Carlo samples.

With 300 exact-model repetitions and the minimum number of bootstrap samples
needed to resolve `alpha_n`, the selection counts were:

| n | bootstrap gamma | asymptotic gamma | bootstrap beta | asymptotic beta |
|---:|---:|---:|---:|---:|
| 75 | 1 | 4 | 4 | 10 |
| 150 | 1 | 4 | 1 | 10 |
| 300 | 0 | 3 | 0 | 12 |

No bootstrap refits failed. Null MISE also improved slightly because the rare
false expansions are expensive, especially for beta. Calibration is therefore
doing what it should.

The hard decision still loses too much power. Two extra scenarios make this
visible: a closer gamma mixture and an interior beta mixture whose fitted beta
shapes remain above the boundary guard. For the stronger mixtures:

| n | method | gamma-mixture selections / MISE | interior-beta selections / MISE |
|---:|---|---:|---:|
| 75 | always use degree 8 | 100 / 0.00781 | 100 / 0.09049 |
|  | asymptotic shrinking gate | 4 / 0.01155 | 99 / 0.09411 |
|  | bootstrap shrinking gate | 3 / 0.01155 | 40 / 0.37528 |
| 150 | always use degree 8 | 100 / 0.00625 | 100 / 0.04450 |
|  | asymptotic shrinking gate | 17 / 0.01015 | 100 / 0.04450 |
|  | bootstrap shrinking gate | 1 / 0.01097 | 86 / 0.11659 |
| 300 | always use degree 8 | 100 / 0.00397 | 100 / 0.02795 |
|  | asymptotic shrinking gate | 79 / 0.00504 | 100 / 0.02795 |
|  | bootstrap shrinking gate | 11 / 0.00966 | 100 / 0.02795 |

The bootstrap methods cost about 0.08--0.40 seconds per fit in these runs,
versus roughly 0.002--0.004 seconds for the analytic gates. Per-fit bootstrap
is therefore useful as a calibration experiment, not as the speed-oriented
end product. Reusable critical-value tables or a response-surface
approximation could retain the finite-sample correction cheaply.

Second, the beta experiment now declines the beta-start route when either
fitted shape is at most 1.1. A cheap endpoint diagnostic based on the first
`floor(sqrt(n))` order statistics then uses uniform-start UCV when an endpoint
power is estimated at most one, and uniform-start HS otherwise. In the paired
50-replication pilot this:

- guarded all 50 `Beta(0.8,0.8)` samples and sent 46 to UCV;
- guarded all 50 beta-mixture samples and sent all 50 to HS;
- did not guard the exact-beta, Kumaraswamy, or logistic-normal samples.

A fresh 500-replication run guarded 499 boundary samples and all 500 mixture
samples. The endpoint diagnostic sent 458 boundary samples and 10 mixture
samples to UCV, using HS for the rest. MISE was 0.03900 on the boundary target
and 0.04653 on the mixture; the latter compares with about 0.31 for the
parametric beta fit. It never activated on the exact-beta, Kumaraswamy, or
logistic-normal samples, whose MISE was 0.01241, 0.01104, and 0.01741. The
conservative global gate still smoothed 11/500 exact-beta samples, however,
and those rare samples account for a disproportionate share of its risk.

Both experiment scripts now use midpoint quadrature in their boundary-
resolving transformed coordinate. For 100 identical `Beta(0.8,0.8)` samples,
the guarded MISE estimates were 0.03651, 0.03648, and 0.03646 with 300, 500,
and 1,000 points. The older endpoint trapezoids were not converged and their
boundary MISE values should not be used. The shape and endpoint cutoffs were
developed on the first small pilot, so the fresh run is encouraging but not
enough to make them public defaults.

## 6. A continuous transition

The hard bootstrap decision solves the null-tail problem but throws away too
much information immediately below its cutoff. The new experimental methods
`gamma_LaguerreK_bootstrap_smooth` and
`beta_JacobiK_bootstrap_smooth_guarded` instead use

```text
w_p(T, c) = (T / c)^p / {1 + (T / c)^p},
```

where `T` is the observed global statistic and `c` is its fitted-null
bootstrap critical value. For gamma, the weight scales the estimated
roughness and therefore enters the analytic bandwidth directly:

```text
J_smooth = w J,
b_smooth = b_Laguerre / w^(2/5).
```

This continuously joins the finite analytic bandwidth to the infinite-
bandwidth parametric endpoint. The beta kernel implementation restricts
finite bandwidths to less than `1/4`, so it cannot follow that path to
infinity. Its smooth method instead averages the parametric density and the
finite-bandwidth Jacobi density with weight `w`. The beta boundary guard
remains hard: it is a safety condition, not a model-selection cutoff.

A transition-power comparison at `n = 150` used 100 repetitions and 149
bootstrap refits per fit. Gamma favored the gentle `p = 1` transition; beta
favored `p = 4`, which better separated the interior mixture from
near-parametric targets. Follow-up pilots used 100 repetitions at `n = 75`
and `150`, and 50 at `n = 300`. Entries below are
parametric / smooth / always-on degree-eight MISE:

| n | exact gamma | gamma mixture | lognormal |
|---:|---:|---:|---:|
| 75 | .00419 / .00564 / .00649 | .01213 / .01095 / .00823 | .01294 / .00908 / .00804 |
| 150 | .00174 / .00225 / .00264 | .01063 / .00962 / .00589 | .01063 / .00767 / .00513 |
| 300 | .00093 / .00127 / .00143 | .01005 / .00700 / .00433 | .01001 / .00344 / .00297 |

| n | exact beta | interior beta mixture | logistic-normal |
|---:|---:|---:|---:|
| 75 | .02417 / .02578 / .05578 | .56006 / .23117 / .08342 | .03154 / .03252 / .06336 |
| 150 | .01127 / .01146 / .02980 | .55522 / .07310 / .04938 | .01665 / .01681 / .03184 |
| 300 | .00504 / .00504 / .01423 | .54719 / .02937 / .02650 | .01112 / .01051 / .01511 |

The continuous rule is not oracle-optimal, but it removes the bootstrap
gate's worst behavior: its interior-beta mixture MISE was 0.375, 0.117, and
0.028 under the hard gate, compared with 0.231, 0.073, and 0.029 here.
For exact beta it is essentially parametric by `n = 300`; for fixed
alternatives it moves toward the finite-bandwidth estimator with sample size.
Gamma shows the same interpolation more gradually. These are tuning pilots,
not final accuracy claims: the transition powers were selected on the
reported scenarios and the `n = 300` run has only 50 repetitions.

## 7. Why beta starts need the separate guard

The beta boundary scenario deserves separate treatment. If a fitted beta
shape is at most `1/2`, the squared fitted density is not integrable at that
endpoint. Even above `1/2`, estimated exponents close to the boundary can make
finite-sample ISE extremely unstable. A finite evaluation grid hides a truly
infinite ISE when a fitted shape crosses `1/2`.

This is not primarily a bandwidth-selection problem. Any eventual public
Jacobi selector should decline the beta-start route in the boundary regime
and use the existing uniform-start HS/UCV machinery. The current research
helper's `shape > 1/2` condition is the mathematical minimum, not a recommended
production threshold.

## 8. Current decision

The work now supports a fairly clear product decision:

1. Keep `GR` and `HS` as the cheap uniform-start reference rules.
2. Keep UCV as the robust public selector for nonuniform asymmetric starts.
3. Do not yet expose a fixed-degree Laguerre or Jacobi rule.
4. Continue with the smooth transition, but replace its per-fit bootstrap by
   pretabulated calibration before considering a public implementation; do not
   use the hard bootstrap gate as a default.
5. Treat the beta boundary regime as an explicit fallback case, not something
   a Jacobi expansion should repair.

The analytic selectors remain valuable: they are fast experimental probes,
and the oracle study confirms that the underlying estimator families can beat
uniform-start smoothing. What is missing is a stable data-driven transition
between the parametric null and complex alternatives.
