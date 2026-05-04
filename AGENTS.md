# AGENTS.md

This repository is an R package. Keep changes small, package-oriented, and easy to verify locally.

## Core Commands

- `just test`: run the testthat suite.
- `just check`: run `R CMD check` with `--as-cran` and `--no-manual`.
- `just coverage`: print coverage and write both `coverage.xml` and a readable `coverage.html`.
- `just readme`: render `README.md` from `README.Rmd`.
- `just clean`: remove local check, coverage, tarball, and generated test artifacts.

## Repo Conventions

- `README.md` is generated from `README.Rmd`. Edit `README.Rmd` if the rendered README needs to change.
- Package documentation in `man/` is generated from roxygen comments in `R/`.
- Keep top-level developer files such as `AGENTS.md` and `Justfile` out of package builds.
- `tests/testthat/Rplots.pdf` is a disposable artifact. Recipes should remove it automatically.

## Change Scope

- Prefer focused edits over wide refactors.
- Add or extend tests for behavior changes.
- Before a version bump, make sure `DESCRIPTION`, docs, and CI reflect the current package behavior.
