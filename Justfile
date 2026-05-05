set shell := ["bash", "-cu"]

default:
    @just --list

test:
    trap 'rm -f tests/testthat/Rplots.pdf' EXIT; Rscript -e "testthat::test_local(reporter = 'summary')"

check:
    trap 'rm -f tests/testthat/Rplots.pdf' EXIT; pkg=$(Rscript -e 'd <- read.dcf("DESCRIPTION"); cat(sprintf("%s_%s.tar.gz", d[1, "Package"], d[1, "Version"]))'); R CMD build . && R CMD check --as-cran --no-manual "$pkg"

coverage:
    trap 'rm -f tests/testthat/Rplots.pdf' EXIT; Rscript tools/render-coverage-report.R

revdep:
    R_COMPILE_AND_INSTALL_PACKAGES=never Rscript -e "revdepcheck::revdep_check(num_workers = 2)"

revdep-reset:
    Rscript -e "revdepcheck::revdep_reset()"

revdep-manual:
    #!/usr/bin/env bash
    set -euo pipefail
    rm -rf revdep-manual
    R CMD INSTALL .
    R_COMPILE_AND_INSTALL_PACKAGES=never Rscript -e 'install.packages(c("RealSurvSim","tscopula"), repos="https://cloud.r-project.org")'
    mkdir -p revdep-manual
    cd revdep-manual
    Rscript -e 'download.packages(c("RealSurvSim","tscopula"), destdir=".", type="source", repos="https://cloud.r-project.org")'
    for tarball in *.tar.gz; do
      echo "=== Checking $tarball ==="
      R CMD check --no-manual "$tarball" || true
    done
    echo
    echo "=== Summary ==="
    grep -H "^Status:" */00check.log

readme:
    Rscript -e "pkgload::load_all(quiet = TRUE); rmarkdown::render('README.Rmd', output_format = 'github_document', quiet = TRUE)"

clean:
    rm -rf ..Rcheck .Rcheck kdensity.Rcheck inst/doc README.html README_files coverage.html coverage.xml cobertura.xml *.tar.gz tests/testthat/Rplots.pdf revdep revdep-manual
