## Update

This release adds compiled evaluation for the built-in kernels, an optional
FFT grid backend, and the explicit `"GR"` gamma-reference bandwidth selector.
It also corrects asymmetric-kernel cross-validation, the vectorized
bias-corrected beta kernel, `logLik()`/AIC/BIC, and `update()` copy semantics.

## Test environments
* local: Ubuntu 26.04, R 4.6.1

## R CMD check results
0 ERRORs, 0 WARNINGs, 0 NOTEs.

## Reverse dependencies
We checked the 2 reverse imports (RealSurvSim and tscopula) by running
R CMD check against this version. No problems were found.
