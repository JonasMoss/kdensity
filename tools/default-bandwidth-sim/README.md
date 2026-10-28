# Default-bandwidth simulation for the asymmetric kernels

Board task: TASK-4 (decisions 2 and 3 of the 2.0.0 plan, 2026-10-28).

The asymmetric kernels default to `bw = "ucv"`, which costs O(n^2) per
optimizer step (48--198 s at n = 10,000). The question is whether the package's
existing closed-form rules can be the defaults instead, applied regardless of
the parametric start (as `RHE` already is for symmetric kernels):

| kernel | rule | candidates in the simulation |
|---|---|---|
| `gamma`, `gamma_biased` | gamma reference (`GR`) | `GR`, `GRc2`, `GRc3` |
| `beta`, `beta_biased` | Hallberg Szabadvary (`HS`) | `HS` |
| `gcopula` | Jones--Henderson (`JH`) | `JH` |

`GR` is undefined when the fitted gamma shape is at most 3/2. `GRc2` and
`GRc3` cap the fitted shape below at 2 and 3 while keeping the fitted mean;
they equal `GR` above the cap. The biased kernels reuse the corrected
kernels' rules.

Each replication also records the package's `ucv` (`ucv_pkg`, current search
interval), UCV minimized over a wide log interval around the rule
(`ucv_wide`, `[rule / 10, 10 rule]` capped at 1/4 for beta kernels and 1/2 for
`gcopula`), and the ISE oracle over the same interval. All runs use the fixed
gamma kernel and the normalized UCV criterion.

Design: 8 kernel x start combinations, 7 targets on the positive half-line or
4 on the unit interval, n = 100 and 500 (200 replications) and n = 2000 (100),
in blocks of 50 replications: 410 cells. Data for a given target, n and
replication are identical across kernels and starts.

## Running

```sh
# Package with the fixed gamma kernel installed into a private library.
R CMD INSTALL --library="$KDENSITY_LIB" <kdensity checkout>
KDENSITY_LIB=... ./run.sh 12      # 12 parallel workers under nice -n 10
Rscript combine.R results          # summary.csv + MISE-ratio table
```

`results/<cell>.csv` is written atomically and marks a finished cell, so
`run.sh` resumes after an interruption. GNU xargs lowers its parallelism by
one per `SIGUSR2` and raises it per `SIGUSR1`.

On simbox (2026-10-28): `~/kdensity-sim/{pkg,lib,sim}`; launched with
16 workers, lowered to 12 because another job shared the machine.

## Results (2026-10-28)

All 410 cells finished without failures. Entries are MISE divided by the
oracle MISE: the geometric mean over targets x n, then the maximum. Selection
time is the mean at n = 2000. `summary.csv` has every cell;
`probit-rules.csv` has the Gaussian-copula follow-up (`probit-rules.R`).

| kernel | start | rule | rule: gmean / max | `ucv` (package): gmean / max | time rule / `ucv` |
|---|---|---|---|---|---|
| `gamma` | uniform | GRc3 | 1.23 / 1.93 | 1.37 / 2.62 | 0.01 s / 1.9 s |
| `gamma` | gamma | GRc3 | 1.67 / 3.25 | 1.38 / 2.11 | 0.01 s / 4.4 s |
| `gamma_biased` | uniform | GRc3 | 1.28 / 2.29 | 1.40 / 2.83 | 0.02 s / 1.9 s |
| `beta` | uniform | HS | 1.13 / 1.79 | 2.58 / 22.4 | 0.001 s / 2.4 s |
| `beta` | beta | HS | 1.85 / 2.84 | 1.96 / 26.3 | 0.001 s / 7.1 s |
| `beta_biased` | uniform | HS | 1.17 / 1.75 | 2.61 / 17.4 | 0.001 s / 1.8 s |
| `gcopula` | uniform | JH | 1.98 / 29.0 | 1.59 / 3.28 | 0.001 s / 1.7 s |
| `gcopula` | uniform | SJ, probit scale | 1.44 / 2.59 | | |
| `gcopula` | beta | JH | 3.16 / 28.9 | 2.05 / 4.43 | 0.001 s / 2.5 s |
| `gcopula` | beta | SJ, probit scale | 2.52 / 6.10 | | |

Findings:

- GR needs a shape floor. Plain GR is undefined in most exponential, gamma
  mixture and lognormal(0, 1) replications (fitted shape at most 3/2). A floor
  of 3 (GRc3) is the most robust variant; a floor of 2 is worse on the mixture
  and the heavy-tailed lognormal.
- HS works for both beta kernels. Package UCV collapses on the beta mixture
  because its search interval `[JH / 4, 1 / 4]` starts at a huge JH value.
- The start-agnostic rules undersmooth when a parametric start fits well
  (gamma start on gamma/exponential targets, beta start on beta targets), where
  UCV exploits the start. With the rules, the start still lowers absolute MISE
  in 15 of 21 gamma-kernel and 9 of 12 beta-kernel cells.
- JH is a normal-reference rule on the probit scale and fails on multimodal
  data (bimodal beta mixture: 4-29x oracle), with or without a guard at its 0.5
  cap. On the probit scale the Gaussian-copula kernel is a Gaussian kernel with
  sd `s = sqrt(1 - rho^2)`, so a Gaussian bandwidth `s` for `qnorm(x)` maps to
  `h = sqrt(1 - sqrt(1 - s^2))`. `stats::bw.SJ` there is robust; JH remains
  better on unimodal and U-shaped targets.
- Caveat: the oracle searches `[rule / 10, 10 rule]`. When a start fits very
  well the true optimum can lie above that (package UCV beat the oracle for
  beta(0.8, 0.8) with a beta start), so ratios for good starts are slightly
  optimistic for UCV.
