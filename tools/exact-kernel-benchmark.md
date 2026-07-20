# Exact sorted-kernel benchmark

This note records the exact evaluation optimization for the compact-support
and Laplace kernels. It complements the approximate FFT grid backend:

- `kdensity_grid(..., method = "fft")` is intended for fast approximate grids;
- callable `kdensity` fits remain exact and now acquire a faster sorted path
  after enough cumulative work.

## Algorithms

For rectangular, triangular, Epanechnikov, biweight, triweight, tricube,
cosine, and optimal-cosine kernels, evaluation first finds the observations
inside `[y - h, y + h]` by binary search. A segment tree answers the resulting
weighted range queries.

Each usable tree node stores moments in coordinates centered within that
node, scaled by `h`. Polynomial kernels need moments through degree two, four,
six, or nine; the cosine kernels store local sine and cosine sums. Nodes wider
than `2h` are not summarized and are recursively split when queried. Thus a
selected node always has bounded local coordinates. This avoids the severe
cancellation that global raw moments would suffer when, for example, data
near `1e8` are smoothed with `h = 0.2`.

The Laplace kernel is simpler. After sorting data and evaluation points, a
forward recurrence accumulates observations to the left and a backward
recurrence accumulates observations to the right. All exponential updates
have non-positive arguments, avoiding overflow.

Sorting and cache construction are lazy. Compact kernels switch after
cumulative direct work reaches

```text
max(50,000, 150 n)
```

data/evaluation pairs. Laplace uses `max(10,000, 10 n)` because it has no
moment-tree construction. Repeated scalar calls count toward the threshold.

## Correctness

Tests compare all nine specialized evaluators with the original compiled
`O(n m)` summation, including:

- named, missing, infinite, duplicated, and unsorted evaluation points;
- parametric-start importance weights;
- bandwidths from `0.02` to `2` for Laplace;
- compact-kernel data centered at both zero and `1e8`.

Maximum discrepancies in the benchmark below were below `2e-14`. In the
large-offset stress test they were below `5e-14`.

## Speed

Run:

```sh
Rscript tools/exact-kernel-benchmark.R \
  --reps=5 --sizes=100,1000,10000 --points=1000
```

The July 2026 development-machine results for 10,000 observations and 1,000
evaluation points were:

| Kernel | Direct (ms) | First sorted (ms) | First speedup | Cached (ms) | Cached speedup |
|---|---:|---:|---:|---:|---:|
| Epanechnikov | 134 | 20 | 6.7x | 1 | 134x |
| Rectangular | 135 | 20 | 6.8x | 1 | 135x |
| Triangular | 138 | 21 | 6.6x | 1 | 138x |
| Biweight | 147 | 20 | 7.4x | 1 | 147x |
| Cosine | 152 | 19 | 8.0x | 1 | 152x |
| Optimal cosine | 148 | 19 | 7.8x | 1 | 148x |
| Triweight | 148 | 20 | 7.4x | 2 | 74x |
| Tricube | 144 | 21 | 6.9x | 4 | 36x |
| Laplace | 167 | 1 | 167x | 1 | 167x |

At 1,000 observations, first-use compact speedups were about 3.5--6.5x and
cached speedups about 6--34x; Laplace reached 56x and 85x. Even at 100
observations and 1,000 evaluation points the threshold-selected path did not
regress, although absolute timings were only one or two milliseconds.
