# Fast-grid verification report

Generated: 2026-10-28 11:31:57 CET

## Decision

**PASS**: accuracy passed and performance passed.

- Promotion accuracy cases: 126
- Promotion accuracy failures: 0
- Characterization-only accuracy cases: 486
- Promotion performance cases: 12
- Promotion performance failures: 0

## Worst promotion accuracy metrics

- Peak-scaled RMSE: 7.562e-05 (gcopula, interior, start=uniform, n=100, h=0.04)
- Peak-scaled maximum error: 0.001094 (gcopula, boundary_heavy, start=beta, n=100, h=0.04)
- Maximum central relative error: 0.009421 (gcopula, boundary_heavy, start=beta, n=100, h=0.04)
- Relative integrated absolute error: 0.0006159 (gcopula, boundary_heavy, start=beta, n=100, h=0.04)
- Relative mass error: 1.625e-05 (gcopula, boundary_heavy, start=beta, n=100, h=0.04)

## Explicit compact-FFT characterization

- Worst peak-scaled RMSE: 0.01097 (rectangular, normal, start=normal, n=100, h=0.04722)
- Worst relative integrated absolute error: 0.02215 (rectangular, skewed, start=uniform, n=1000, h=0.03625)

These cases remain explicit opt-ins; automatic plotting uses their exact
sorted evaluators.

## Warmed evaluation timings

| Kernel | Start | n | Exact ms | FFT ms | Speedup | Accepted |
|---|---|---:|---:|---:|---:|:---:|
| gaussian | uniform | 100 | 1.375 | 1.047 | 1.31 | yes |
| epanechnikov | uniform | 100 | 0.2656 | 0.9141 | 0.291 | yes |
| rectangular | uniform | 100 | 0.1973 | 0.9844 | 0.2 | yes |
| triangular | uniform | 100 | 0.4922 | 1.273 | 0.387 | yes |
| biweight | uniform | 100 | 0.5156 | 1.07 | 0.482 | yes |
| cosine | uniform | 100 | 0.3125 | 1.102 | 0.284 | yes |
| optcosine | uniform | 100 | 0.3125 | 1.109 | 0.282 | yes |
| triweight | uniform | 100 | 0.9688 | 1.852 | 0.523 | yes |
| tricube | uniform | 100 | 1.328 | 2.219 | 0.599 | yes |
| laplace | uniform | 100 | 0.2793 | 2.156 | 0.13 | yes |
| gaussian | uniform | 1000 | 20 | 2.234 | 8.95 | yes |
| epanechnikov | uniform | 1000 | 0.75 | 2.188 | 0.343 | yes |
| rectangular | uniform | 1000 | 0.4883 | 2.156 | 0.226 | yes |
| triangular | uniform | 1000 | 0.8125 | 1.656 | 0.491 | yes |
| biweight | uniform | 1000 | 0.832 | 1.641 | 0.507 | yes |
| cosine | uniform | 1000 | 0.5352 | 1.703 | 0.314 | yes |
| optcosine | uniform | 1000 | 0.7656 | 2.172 | 0.353 | yes |
| triweight | uniform | 1000 | 1.391 | 2.266 | 0.614 | yes |
| tricube | uniform | 1000 | 2.594 | 2.313 | 1.12 | yes |
| laplace | uniform | 1000 | 0.4141 | 1.875 | 0.221 | yes |
| gaussian | uniform | 10000 | 202 | 3.781 | 53.4 | yes |
| epanechnikov | uniform | 10000 | 1.234 | 3.688 | 0.335 | yes |
| rectangular | uniform | 10000 | 0.8828 | 6.75 | 0.131 | yes |
| triangular | uniform | 10000 | 1.758 | 3.344 | 0.526 | yes |
| biweight | uniform | 10000 | 1.656 | 4.031 | 0.411 | yes |
| cosine | uniform | 10000 | 1.094 | 3.625 | 0.302 | yes |
| optcosine | uniform | 10000 | 1.148 | 3.656 | 0.314 | yes |
| triweight | uniform | 10000 | 1.703 | 2.859 | 0.596 | yes |
| tricube | uniform | 10000 | 3.219 | 7.688 | 0.419 | yes |
| laplace | uniform | 10000 | 0.7734 | 2.25 | 0.344 | yes |
| gaussian | normal | 100 | 1.508 | 1.195 | 1.26 | yes |
| gaussian | normal | 1000 | 14.5 | 1.781 | 8.14 | yes |
| gaussian | normal | 10000 | 208 | 3.938 | 52.8 | yes |
| gcopula | uniform | 100 | 4.094 | 3.781 | 1.08 | yes |
| gcopula | beta | 100 | 4.062 | 3.391 | 1.2 | yes |
| gcopula | uniform | 1000 | 37.25 | 14.12 | 2.64 | yes |
| gcopula | beta | 1000 | 38.75 | 14.62 | 2.65 | yes |
| gcopula | uniform | 10000 | 379 | 125 | 3.03 | yes |
| gcopula | beta | 10000 | 451 | 128 | 3.52 | yes |

Full case-level results are in `accuracy.csv` and `performance.csv`.
