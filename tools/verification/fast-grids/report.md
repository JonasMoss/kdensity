# Fast-grid verification report

Generated: 2026-10-28 09:06:26 CET

## Decision

**PASS**: accuracy passed and performance passed.

- Promotion accuracy cases: 126
- Promotion accuracy failures: 0
- Characterization-only accuracy cases: 324
- Promotion performance cases: 12
- Promotion performance failures: 0

## Worst promotion accuracy metrics

- Peak-scaled RMSE: 7.783e-05 (gcopula, interior, start=uniform, n=100, h=0.04)
- Peak-scaled maximum error: 0.001127 (gcopula, boundary_heavy, start=beta, n=100, h=0.04)
- Maximum central relative error: 0.009859 (gcopula, boundary_heavy, start=beta, n=100, h=0.04)
- Relative integrated absolute error: 0.0006373 (gcopula, boundary_heavy, start=beta, n=100, h=0.04)
- Relative mass error: 1.851e-05 (gcopula, boundary_heavy, start=beta, n=100, h=0.04)

## Explicit compact-FFT characterization

- Worst peak-scaled RMSE: 0.01511 (rectangular, normal, start=normal, n=100, h=0.04722)
- Worst relative integrated absolute error: 0.02052 (rectangular, normal, start=uniform, n=5000, h=0.05037)

These cases remain explicit opt-ins; automatic plotting uses their exact
sorted evaluators.

## Warmed evaluation timings

| Kernel | Start | n | Exact ms | FFT ms | Speedup | Accepted |
|---|---|---:|---:|---:|---:|:---:|
| gaussian | uniform | 100 | 2.016 | 2 | 1.01 | yes |
| epanechnikov | uniform | 100 | 0.4336 | 1.828 | 0.237 | yes |
| rectangular | uniform | 100 | 0.3027 | 1.781 | 0.17 | yes |
| triangular | uniform | 100 | 0.5312 | 1.891 | 0.281 | yes |
| biweight | uniform | 100 | 0.6445 | 1.953 | 0.33 | yes |
| cosine | uniform | 100 | 0.4258 | 1.828 | 0.233 | yes |
| optcosine | uniform | 100 | 0.4062 | 1.937 | 0.21 | yes |
| gaussian | uniform | 1000 | 17 | 1.469 | 11.6 | yes |
| epanechnikov | uniform | 1000 | 0.5703 | 1.445 | 0.395 | yes |
| rectangular | uniform | 1000 | 0.3984 | 1.383 | 0.288 | yes |
| triangular | uniform | 1000 | 0.7969 | 1.469 | 0.543 | yes |
| biweight | uniform | 1000 | 0.8437 | 1.414 | 0.597 | yes |
| cosine | uniform | 1000 | 0.5703 | 2.078 | 0.274 | yes |
| optcosine | uniform | 1000 | 0.543 | 1.477 | 0.368 | yes |
| gaussian | uniform | 10000 | 155 | 1.641 | 94.5 | yes |
| epanechnikov | uniform | 10000 | 0.9609 | 1.672 | 0.575 | yes |
| rectangular | uniform | 10000 | 0.7383 | 1.625 | 0.454 | yes |
| triangular | uniform | 10000 | 1.258 | 1.688 | 0.745 | yes |
| biweight | uniform | 10000 | 1.273 | 1.641 | 0.776 | yes |
| cosine | uniform | 10000 | 0.8828 | 1.688 | 0.523 | yes |
| optcosine | uniform | 10000 | 0.9063 | 1.781 | 0.509 | yes |
| gaussian | normal | 100 | 1.563 | 1.477 | 1.06 | yes |
| gaussian | normal | 1000 | 15.5 | 1.453 | 10.7 | yes |
| gaussian | normal | 10000 | 155 | 1.641 | 94.5 | yes |
| gcopula | uniform | 100 | 2.703 | 2.484 | 1.09 | yes |
| gcopula | beta | 100 | 2.875 | 2.656 | 1.08 | yes |
| gcopula | uniform | 1000 | 26.75 | 10 | 2.68 | yes |
| gcopula | beta | 1000 | 27.25 | 9.875 | 2.76 | yes |
| gcopula | uniform | 10000 | 286 | 86.5 | 3.31 | yes |
| gcopula | beta | 10000 | 281 | 83 | 3.39 | yes |

Full case-level results are in `accuracy.csv` and `performance.csv`.
