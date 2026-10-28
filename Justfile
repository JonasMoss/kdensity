set shell := ["bash", "-cu"]

default:
    @just --list

test:
    trap 'rm -f tests/testthat/Rplots.pdf' EXIT; Rscript -e "testthat::test_local(reporter = 'summary')"

check:
    trap 'rm -f tests/testthat/Rplots.pdf' EXIT; pkg=$(Rscript -e 'd <- read.dcf("DESCRIPTION"); cat(sprintf("%s_%s.tar.gz", d[1, "Package"], d[1, "Version"]))'); R CMD build . && R CMD check --as-cran --no-manual "$pkg"

coverage:
    trap 'rm -f tests/testthat/Rplots.pdf' EXIT; Rscript tools/render-coverage-report.R

verify-grids:
    Rscript tools/verification/fast-grids/verify.R

revdep:
    R_COMPILE_AND_INSTALL_PACKAGES=never Rscript -e "revdepcheck::revdep_check(num_workers = 2)"

revdep-reset:
    Rscript -e "revdepcheck::revdep_reset()"

revdep-manual:
    #!/usr/bin/env bash
    set -euo pipefail
    # Private library with Posit Package Manager binaries for this Linux release;
    # the revdeps' Suggests are installed too, since R CMD check needs them.
    revdeps='c("RealSurvSim", "tscopula", "TreeDist")'
    . /etc/os-release
    p3m="https://packagemanager.posit.co/cran/__linux__/${VERSION_CODENAME}/latest"
    rm -rf revdep-manual
    mkdir -p revdep-manual/lib
    export R_LIBS_USER="$PWD/revdep-manual/lib"
    Rscript -e "install.packages($revdeps, dependencies = TRUE, repos = '$p3m', Ncpus = 4)"
    R CMD INSTALL .
    cd revdep-manual
    Rscript -e "download.packages($revdeps, destdir = '.', type = 'source', repos = 'https://cloud.r-project.org')"
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
