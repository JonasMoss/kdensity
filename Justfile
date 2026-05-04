set shell := ["bash", "-cu"]

default:
    @just --list

test:
    Rscript -e "testthat::test_local(reporter = 'summary')"

check:
    pkg=$(Rscript -e 'd <- read.dcf("DESCRIPTION"); cat(sprintf("%s_%s.tar.gz", d[1, "Package"], d[1, "Version"]))'); R CMD build . && R CMD check --as-cran --no-manual "$pkg"

coverage:
    Rscript -e "cov <- covr::package_coverage(); print(cov); cat(as.character(covr::to_cobertura(cov)), file = 'coverage.xml'); if (requireNamespace('DT', quietly = TRUE) && requireNamespace('htmltools', quietly = TRUE)) { covr::report(cov, file = 'coverage.html', browse = FALSE) } else { message('Skipping coverage.html; install DT and htmltools to render it.') }"
