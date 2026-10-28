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
