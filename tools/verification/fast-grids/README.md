# Fast-grid verification

This directory contains promotion-grade verification for the approximate FFT
grid backend. The entire `tools/` tree is excluded from package builds by
`.Rbuildignore`.

Run the full verification from the repository root:

```sh
just verify-grids
```

For a quicker development pass:

```sh
Rscript tools/verification/fast-grids/verify.R --quick
```

The script writes:

- `accuracy.csv`: one row per accuracy case;
- `performance.csv`: warmed exact-versus-FFT timings;
- `report.md`: acceptance summary and promotion decision.

## Accuracy contract

The matrix characterizes every explicit FFT path and applies the promotion
contract to the automatic Gaussian and Gaussian-copula routes. Compact kernels
remain available through explicit `method = "fft"`, but `"auto"` keeps their
faster and more accurate exact sorted evaluators.

The matrix covers:

- every supported translation-invariant kernel;
- the Gaussian-copula kernel;
- uniform and representative parametric starts;
- normal, skewed, multimodal, interior-beta, and boundary-heavy beta targets;
- small, medium, and large samples;
- narrow, ordinary, and wide bandwidths;
- the default plotting grid size and default internal FFT resolution.

The broad compact-kernel/parametric-start matrix uses `normalized = FALSE` to
isolate the weighted convolution from the unrelated numerical-integration
limitations of discontinuous kernels. Normalized Gaussian/normal and
Gaussian-copula/beta paths remain covered by the package tests and performance
cases; multiplying both exact and FFT results by the same fitted normalization
does not change their relative approximation error.

Promotion cases must return finite, non-negative values wherever exact
evaluation is finite and satisfy:

- peak-scaled RMSE no larger than 0.2%;
- peak-scaled maximum error no larger than 1%;
- maximum relative error where the exact density exceeds 1% of its peak no
  larger than 2%;
- relative integrated absolute error and relative mass error no larger than
  0.3%.

Explicit compact-kernel FFT cases are retained as characterization data. Their
relative error near compact-support zeros, and uniform pointwise convergence
for the discontinuous rectangular kernel, are not used to promote the
automatic route.

## Performance contract

Performance is measured after warming both evaluators. Promotion requires:

- at least a two-fold speedup for Gaussian and Gaussian-copula grids at
  `n >= 1,000`;
- small-sample FFT evaluation to remain below 10 ms even when its ratio is
  dominated by timer and dispatch overhead.

Compact-kernel timings verify that `"auto"` should retain exact routing when
the sorted evaluator is faster than FFT.
