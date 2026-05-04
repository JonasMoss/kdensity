## Resubmission / update

* New closed-form bandwidth selector `"HS"` (Hallberg Szabadváry's beta
  reference rule, arXiv:2601.19553); now the default for the
  bias-corrected `beta` kernel with a uniform or constant start.
* Several packages removed from Imports and Suggests.
* Internal cleanup: consolidated R source layout, no user-visible API
  changes.

## Test environments
* macos-latest (release)
* windows-latest (release)
* ubuntu-latest (devel)
* ubuntu-latest (release)
* ubuntu-latest (oldrel-1)

## R CMD check results
0 ERRORs, 0 WARNINGs, 2 NOTEs:

* `Found the following (possibly) invalid URLs:` — the Biometrika article
  link in `README.md` returns HTTP 403 to automated requests but is
  reachable in a browser. The DOI is correct.
* `unable to verify current time` — local clock-verification service
  unavailable; not package-related.

## Reverse dependencies
No problems.
