context("utility functions")

expect_silent(assert_(TRUE))
expect_silent(assert_(TRUE, TRUE))
expect_error(assert_(FALSE), "Assertion failed.")
expect_error(assert_(TRUE, FALSE, msg = "custom message"), "custom message")
