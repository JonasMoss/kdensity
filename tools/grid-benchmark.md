# FFT grid benchmark

This note records the focused grid-evaluation optimization implemented in
`R/grid.R`. Callable `kdensity` fits still use exact compiled summation. The
new `kdensity_grid()` function, and `grid = "fft"` in the plotting methods,
explicitly opt into linear binning followed by an FFT convolution.

## Supported kernels

The Gaussian, rectangular, triangular, Epanechnikov, biweight, cosine, and
optimal-cosine kernels are translation invariant and use the convolution
implemented by `stats::density()`. The parametric-start correction is retained
by binning masses proportional to the reciprocal fitted start density.

The Gaussian-copula kernel has a less obvious convolution. With

```text
z = qnorm(y),
q_i = qnorm(x_i),
rho = 1 - h^2,
sigma^2 = 1 - rho^2,
```

its kernel can be written exactly as

```text
K(y, x_i)
  = exp(z^2 / 2) sqrt(2 pi) sigma
    phi_sigma(z - rho q_i).
```

Thus the correction is an ordinary weighted Gaussian convolution on a probit
grid. Literal endpoints zero and one are evaluated exactly; the FFT covers
the finite padded probit interval between them.

Gamma and beta kernels remain on the exact compiled path. Their shape
parameters change with the evaluation point, so an ordinary convolution does
not apply.

## Accuracy and speed

Run the benchmark with:

```sh
Rscript tools/grid-benchmark.R \
  --reps=5 --sizes=100,1000,10000 --grid_n=1000 --fft_n=4096
```

The July 2026 development-machine run gave the following median evaluation
times. Construction and bandwidth selection are excluded.

| Case | Sample size | Exact (ms) | FFT (ms) | Speedup | Maximum absolute error |
|---|---:|---:|---:|---:|---:|
| Gaussian, uniform start | 100 | 1 | 1 | 1x | 6.5e-6 |
| Gaussian, normal start | 100 | 1 | 1 | 1x | 7.0e-6 |
| Gaussian copula, uniform start | 100 | 2 | 2 | 1x | 3.6e-5 |
| Gaussian copula, beta start | 100 | 2 | 1 | 2x | 3.7e-5 |
| Gaussian, uniform start | 1,000 | 14 | 1 | 14x | 2.3e-6 |
| Gaussian, normal start | 1,000 | 14 | 2 | 7x | 2.1e-6 |
| Gaussian copula, uniform start | 1,000 | 22 | 2 | 11x | 1.6e-5 |
| Gaussian copula, beta start | 1,000 | 21 | 1 | 21x | 7.8e-6 |
| Gaussian, uniform start | 10,000 | 141 | 1 | 141x | 1.1e-6 |
| Gaussian, normal start | 10,000 | 143 | 1 | 143x | 3.7e-7 |
| Gaussian copula, uniform start | 10,000 | 219 | 2 | 110x | 1.3e-5 |
| Gaussian copula, beta start | 10,000 | 221 | 2 | 111x | 7.3e-6 |

At `n = 100` there is no material timing gain, but there is also little
overhead. The benefit becomes large for samples of 1,000 or more. Errors are
controlled by `fft_n` and Gaussian tail padding by `truncate`. Discontinuous
rectangular kernels have larger pointwise error close to their jumps, so their
tests use integrated root-mean-square error rather than an unrealistically
strict uniform bound.
