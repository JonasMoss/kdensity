context("utility functions")

a <- 1:3
b <- letters[2:9]
c <- 9:20
set.seed(313)
expect_equal(length(recycle(a, b, c, prototype = "c")$b), length(c))
expect_equal(length(recycle(a, b, c, prototype = c)$a), length(c))
expect_equal(length(recycle(a, b, c, prototype = 5)$c), 5)
expect_equal(length(recycle(a, b, c, prototype = rnorm(100))$c), 100)
expect_equal(length(recycle(a, b, c)$c), max(length(a), length(b), length(c)))
expect_error(recycle(a, b, c, prototype = function() NULL))
expect_error(recycle(a, b, c, prototype = "haha"))
expect_error(recycle(a, b, c, prototype = -1))

x <- list(
  a = 5,
  b = 0,
  c = "a",
  d = NULL
)

y <- list(
  a = 3,
  b = 7,
  f = NA
)

expect_equal(listmerge(x, y, type = "merge")$b, 7)
expect_null(listmerge(x, y, type = "template")$f, NULL)
expect_equal(listmerge(x, NULL, type = "merge"), x)
expect_equal(listmerge(x, list(h = "b"), type = "merge")$h, "b")
