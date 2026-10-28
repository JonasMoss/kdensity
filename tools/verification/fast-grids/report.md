# Fast-grid verification report

Generated: 2026-10-28 12:32:23 CET

## Decision

**PASS**: accuracy passed and performance passed.

- Promotion accuracy cases: 126
- Promotion accuracy failures: 0
- Characterization-only accuracy cases: 486
- Promotion performance cases: 12
- Promotion performance failures: 0

## Worst promotion accuracy metrics

- Peak-scaled RMSE: 0.000104 (gcopula, boundary_heavy, start=uniform, n=5000, h=0.04)
- Peak-scaled maximum error: 0.001302 (gcopula, boundary_heavy, start=uniform, n=1000, h=0.04)
- Maximum central relative error: 0.01093 (gcopula, boundary_heavy, start=beta, n=100, h=0.04)
- Relative integrated absolute error: 0.000653 (gcopula, boundary_heavy, start=beta, n=100, h=0.04)
- Relative mass error: 1.867e-05 (gcopula, boundary_heavy, start=beta, n=100, h=0.04)

## Explicit compact-FFT characterization

- Worst peak-scaled RMSE: 0.01097 (rectangular, normal, start=normal, n=100, h=0.04722)
- Worst relative integrated absolute error: 0.02215 (rectangular, skewed, start=uniform, n=1000, h=0.03625)

These cases remain explicit opt-ins; automatic plotting uses their exact
sorted evaluators.

## Warmed evaluation timings

| Kernel | Start | n | Exact ms | FFT ms | Speedup | Accepted |
|---|---|---:|---:|---:|---:|:---:|
| gaussian | uniform | 100 | 2.828 | 2.047 | 1.38 | yes |
| epanechnikov | uniform | 100 | 0.4922 | 1.625 | 0.303 | yes |
| rectangular | uniform | 100 | 0.3535 | 1.703 | 0.208 | yes |
| triangular | uniform | 100 | 0.6562 | 1.641 | 0.4 | yes |
| biweight | uniform | 100 | 0.7109 | 1.656 | 0.429 | yes |
| cosine | uniform | 100 | 0.5352 | 1.672 | 0.32 | yes |
| optcosine | uniform | 100 | 0.4648 | 2.313 | 0.201 | yes |
| triweight | uniform | 100 | 1.047 | 1.977 | 0.53 | yes |
| tricube | uniform | 100 | 1.242 | 2.422 | 0.513 | yes |
| laplace | uniform | 100 | 0.4023 | 1.812 | 0.222 | yes |
| gaussian | uniform | 1000 | 27.25 | 1.578 | 17.3 | yes |
| epanechnikov | uniform | 1000 | 0.8125 | 1.75 | 0.464 | yes |
| rectangular | uniform | 1000 | 0.7031 | 1.859 | 0.378 | yes |
| triangular | uniform | 1000 | 1.234 | 1.734 | 0.712 | yes |
| biweight | uniform | 1000 | 1.344 | 1.672 | 0.804 | yes |
| cosine | uniform | 1000 | 0.9141 | 1.906 | 0.48 | yes |
| optcosine | uniform | 1000 | 0.8203 | 1.812 | 0.453 | yes |
| triweight | uniform | 1000 | 1.703 | 1.969 | 0.865 | yes |
| tricube | uniform | 1000 | 3.109 | 2.25 | 1.38 | yes |
| laplace | uniform | 1000 | 0.4219 | 1.766 | 0.239 | yes |
| gaussian | uniform | 10000 | 269 | 2.312 | 116 | yes |
| epanechnikov | uniform | 10000 | 1.547 | 2.234 | 0.692 | yes |
| rectangular | uniform | 10000 | 1.055 | 1.945 | 0.542 | yes |
| triangular | uniform | 10000 | 1.469 | 1.195 | 1.23 | yes |
| biweight | uniform | 10000 | 1.484 | 1.375 | 1.08 | yes |
| cosine | uniform | 10000 | 1 | 1.32 | 0.757 | yes |
| optcosine | uniform | 10000 | 1.023 | 1.32 | 0.775 | yes |
| triweight | uniform | 10000 | 2.062 | 1.328 | 1.55 | yes |
| tricube | uniform | 10000 | 4.094 | 1.445 | 2.83 | yes |
| laplace | uniform | 10000 | 1 | 1.211 | 0.826 | yes |
| gaussian | normal | 100 | 1.812 | 1.25 | 1.45 | yes |
| gaussian | normal | 1000 | 17.88 | 1.305 | 13.7 | yes |
| gaussian | normal | 10000 | 179 | 1.359 | 132 | yes |
| gcopula | uniform | 100 | 3.437 | 1.437 | 2.39 | yes |
| gcopula | beta | 100 | 3.719 | 1.664 | 2.23 | yes |
| gcopula | uniform | 1000 | 67 | 2.328 | 28.8 | yes |
| gcopula | beta | 1000 | 84.5 | 2.547 | 33.2 | yes |
| gcopula | uniform | 10000 | 629 | 4.625 | 136 | yes |
| gcopula | beta | 10000 | 838 | 3.563 | 235 | yes |

Full case-level results are in `accuracy.csv` and `performance.csv`.
