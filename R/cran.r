# To avoid those pesky CRAN notes.

g <- function() {
  eql <- EQL::hermite(1, 1)
  uml <- univariateML::egypt
  c(eql, uml)
}
