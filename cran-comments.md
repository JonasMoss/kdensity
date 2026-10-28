## Submission

This is a major release with breaking changes, listed in NEWS.md: the gamma
kernel now follows Chen (2000), every built-in kernel defaults to a
closed-form bandwidth rule, and fitted objects are read-only. The reverse
dependencies are not affected.

It is the first CRAN version of kdensity with compiled code (C++17, via
Rcpp). We checked the compiled code with AddressSanitizer, UBSan and
valgrind.

It also fixes the r-devel NOTE on the current CRAN version about
'bandwidths.Rd' documenting arguments without \usage.

## Test environments

* local: Ubuntu 26.04, R 4.6.1
* local, sanitizers: gcc 15 AddressSanitizer and UBSan (Ubuntu 26.04,
  R 4.6.1); valgrind (Ubuntu 24.04, R 4.6.1)
* win-builder: R-devel, R-release
* macOS builder: R-release (arm64)
* GitHub Actions: macOS and Windows (R-release); Ubuntu (R-devel, R-release,
  R-oldrel-1)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

We ran R CMD check on all three reverse dependencies against this version:
RealSurvSim and tscopula (Imports) and TreeDist (Suggests). No problems were
found.
