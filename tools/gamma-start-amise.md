# Gamma-kernel AMISE with a parametric start

This note derives an experimental non-cross-validated bandwidth for the
bias-corrected gamma kernel used by `kdensity`. It is research code, not yet a
package contract. The starting expansion is from Song Xi Chen,
["Probability Density Function Estimation Using Gamma Kernels"
(2000)](https://doi.org/10.1023/A:1004165218295).

## 1. Leading bias and variance

Let `g` be a fixed parametric start and write the target density as

\[
f(x)=g(x)r(x).
\]

Before final normalization, the population version of the estimator is

\[
\widetilde f_b(x)
=g(x)\int K_{\rho_b(x),b}(u)r(u)\,du.
\]

For \(x\geq 2b\), the bias-corrected gamma kernel has mean \(x\) and variance
\(xb\). Chen's interior expansion therefore gives

\[
E\{\widetilde f_b(x)\}
=f(x)+\frac b2 xg(x)r''(x)+o(b).
\]

Define

\[
D_g(x)=xg(x)r''(x),\qquad A_g=\int_0^\infty D_g(x)\,dx.
\]

Since `kdensity` divides by the integral of the estimated curve, expanding
that normalizer gives

\[
\operatorname{Bias}\{\widehat f_b(x)\}
=\frac b2\{D_g(x)-A_gf(x)\}+o(b).
\]

The start cancels from the leading local variance, leaving Chen's integrated
variance:

\[
\int\operatorname{Var}\{\widehat f_b(x)\}\,dx
=\frac{V_f}{2\sqrt{\pi}\,n\sqrt b}
o\{(n\sqrt b)^{-1}\},
\quad
V_f=\int_0^\infty x^{-1/2}f(x)\,dx.
\]

Consequently, with

\[
J_g=\int_0^\infty\{D_g(x)-A_gf(x)\}^2\,dx,
\]

the working AMISE and its minimizer are

\[
\operatorname{AMISE}(b)
=\frac{b^2}{4}J_g+\frac{V_f}{2\sqrt{\pi}\,n\sqrt b},
\qquad
b_*=\left\{\frac{V_f}{2\sqrt{\pi}\,nJ_g}\right\}^{2/5}.
\]

For a uniform start, \(D_g=xf''\). Integration by parts gives
\(\int xf''=f(0)\), so final normalization has no leading effect whenever
\(f(0)=0\). This includes the gamma references for which Chen's roughness is
finite.

## 2. Why a fitted gamma reference alone is insufficient

If the gamma start is exactly correct, then \(r=1\), hence \(r''=0\) and
\(J_g=0\). The first-order optimal bandwidth is infinite: the estimator should
remain parametric when there is no detectable departure from the start.

A finite adaptive rule must therefore estimate departures of \(r\) from one,
not substitute \(f=g\) into the ordinary gamma-kernel formula.

## 3. Generalized Laguerre expansion

Fit \(g=\operatorname{Gamma}(\alpha,\lambda)\), using the rate
parameterization, and set \(t=\lambda x\). Under \(g\), \(t\) has weight

\[
w_\alpha(t)=\frac{t^{\alpha-1}e^{-t}}{\Gamma(\alpha)}.
\]

Use the orthonormal basis

\[
\phi_k(t)=
\frac{L_k^{(\alpha-1)}(t)}
{\sqrt{\Gamma(k+\alpha)/(k!\Gamma(\alpha))}},
\qquad
r(x)\approx R_K(t)=\sum_{k=0}^K\theta_k\phi_k(t),
\]

where

\[
\widehat\theta_k=\frac1n\sum_{i=1}^n\phi_k(\lambda X_i).
\]

The gamma MLE makes the degree-one coefficient essentially zero. For
\(k\geq2\),

\[
\frac{d^2}{dt^2}L_k^{(\alpha-1)}(t)
=L_{k-2}^{(\alpha+1)}(t),
\]

so the second derivative of the truncated ratio is available exactly.

Writing this derivative polynomial as \(S_K(t)\), we obtain

\[
D_g(x)
=\frac{\lambda^2}{\Gamma(\alpha)}
e^{-t}t^\alpha S_K(t).
\]

Both \(A_g\) and \(J_g\) are then finite sums of gamma integrals. The
implementation in `gamma-start-amise.R` evaluates those polynomial integrals
analytically and estimates \(V_f\) by
\(n^{-1}\sum_i X_i^{-1/2}\).

Because \(\widehat\theta_k\) fluctuates at approximately the
\(n^{-1/2}\) scale even under the parametric model, the experiment also
includes a provisional hard-thresholded variant. It replaces coefficients of
degree two and above by zero when

\[
|\widehat\theta_k| < n^{-1/2}.
\]

If every curvature coefficient is removed, \(J_g=0\), the selected bandwidth
is infinite, and `kdensity` returns the fitted gamma start. This threshold is
a diagnostic arm, not a theoretically calibrated final rule.

## 4. Questions for the simulation

The experiment should establish:

1. which truncation degree gives a useful bias-variance tradeoff;
2. whether coefficient shrinkage is needed under the exact gamma model;
3. sensitivity to fitted shapes near \(1/2\) and to observations near zero;
4. whether the final-normalization term materially improves MISE;
5. performance relative to UCV, Chen's gamma-reference rule, and the fully
   parametric fit under gamma, mixture, lognormal, Weibull, and boundary-heavy
   targets.

The derivation still needs a uniform treatment of the \(x<2b\) boundary
region. Parameter estimation is negligible for the fixed-alternative
roughness at first order, while its near-null effect on coefficient covariance
is treated in `series-selector-findings.md`.

## 5. Initial pilot

A small paired Monte Carlo run used 20 repetitions, \(n=150\), 500 integration
points, and degree \(K=4\). These numbers are diagnostics, not publication
quality results.

| Target | Best or representative methods | Estimated MISE |
|---|---|---:|
| Gamma(shape=3, rate=1.5) | Parametric gamma | 0.00208 |
|  | Thresholded Laguerre-4 | 0.00237 |
|  | Raw / normalized Laguerre-4 | 0.00248 / 0.00257 |
|  | Gamma-start UCV | 0.00431 |
| Gamma mixture | Gamma-start UCV | 0.00364 |
|  | Uniform UCV | 0.00557 |
|  | Raw / normalized Laguerre-4 | 0.00774 / 0.00921 |
|  | Parametric gamma | 0.01089 |
| Lognormal | Raw / normalized Laguerre-4 | 0.00562 / 0.00566 |
|  | Uniform UCV | 0.00675 |
|  | Parametric gamma | 0.01002 |
|  | Gamma-start UCV | 0.02338 |
| Weibull | Parametric gamma | 0.00285 |
|  | Raw / normalized Laguerre-4 | 0.00347 / 0.00349 |
|  | Gamma-start UCV | 0.00359 |

The strict gamma-reference rule failed in all 20 boundary-gamma samples and
in 10 of 20 mixture samples because the fitted shape did not exceed \(3/2\).
Its package-level fallback therefore needs to be reported separately from the
strict analytic rule.

Laguerre selection and fitting generally took about 1--2 ms in this run;
gamma-start UCV took roughly 50--100 ms. The pilot suggests:

- thresholding is useful near the exact parametric model;
- a fixed low-order expansion can miss multimodal departures;
- normalization is theoretically required but did not improve this finite
  truncation in the pilot;
- degree selection or coefficient regularization should be investigated
  before adding a package selector.

## 6. Second-stage study

`series-selector-findings.md` records the shared gamma/beta continuation. It
adds fixed oracle-MISE grids, hard/soft/garrote shrinkage, the analytic
fitted-score covariance correction, global smooth tests, and AIC/BIC degree
selection.

The larger degree pilot confirms that degree 6--8 closes much of the oracle
gap for gamma mixtures and lognormal targets, but no automatic rule dominates:
rules that protect the exact gamma model tend to miss mixture departures at
`n = 150`.

`gamma_LaguerreK_conservative` makes that tradeoff explicit by using the
shrinking global-test level `min(0.01, 1 / n)`. In a fresh 500-replication
degree-eight run at `n = 150`, it selected 13 exact-gamma, 103 gamma-mixture,
208 lognormal, and 3 Weibull samples. Its MISE was 0.00174, 0.00963, 0.00759,
and 0.00334, respectively. The 2.6% exact-model selection rate is larger than
the nominal 0.67%, so the asymptotic chi-squared tail itself needs calibration;
lowering the nominal level alone is not a complete fix. The mixture result is
also worse than the 0.00635 obtained when degree eight was always used in the
earlier paired pilot. The experimental selector should therefore remain
outside the package API.

The parametric-bootstrap version fixes most of the null-tail problem: across
300 exact-gamma repetitions it selected 1, 1, and 0 samples at `n = 75, 150,
300`, compared with 4, 4, and 3 for the asymptotic shrinking gate. It is too
conservative as a hard selector, however. On the gamma mixture it selected
only 3%, 1%, and 11% at those sample sizes; at `n = 300` its MISE was 0.00966
versus 0.00397 for always using degree eight. Bootstrap calibration should be
used to build a smoother or tabulated rule, not exposed directly.
