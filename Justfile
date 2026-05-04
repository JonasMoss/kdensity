set shell := ["bash", "-cu"]

default:
    @just --list

test:
    trap 'rm -f tests/testthat/Rplots.pdf' EXIT; Rscript -e "testthat::test_local(reporter = 'summary')"

check:
    trap 'rm -f tests/testthat/Rplots.pdf' EXIT; pkg=$(Rscript -e 'd <- read.dcf("DESCRIPTION"); cat(sprintf("%s_%s.tar.gz", d[1, "Package"], d[1, "Version"]))'); R CMD build . && R CMD check --as-cran --no-manual "$pkg"

coverage:
    trap 'rm -f tests/testthat/Rplots.pdf' EXIT; Rscript tools/render-coverage-report.R

readme:
    Rscript -e "pkgload::load_all(quiet = TRUE); rmarkdown::render('README.Rmd', output_format = 'github_document', quiet = TRUE)"

clean:
    rm -rf ..Rcheck .Rcheck kdensity.Rcheck inst/doc README.html README_files coverage.html coverage.xml cobertura.xml *.tar.gz tests/testthat/Rplots.pdf
