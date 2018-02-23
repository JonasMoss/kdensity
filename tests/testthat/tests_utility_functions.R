context("utility functions")

a = 1:3
b = letters[2:9]
c = 9:20

expect_equal(length(recycle(a, b, c, prototype = "c")$b), length(c))
expect_equal(length(recycle(a, b, c, prototype = c)$a), length(c))
expect_equal(length(recycle(a, b, c, prototype = 5)$c), 5)

x = list(a = 5,
         b = 0,
         c = "a",
         d = NULL)

y = list(a = 3,
         b = 7,
         f = NA)

expect_equal(listmerge(x, y, type = "merge")$b, 7)
expect_null(listmerge(x, y, type = "template")$f, NULL)
