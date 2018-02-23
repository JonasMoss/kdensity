context("bandwidths")

t = function(x) {}
good = function(x, kernel, start, support) {}

expect_error(add_bw("t", t))
expect_silent(add_bw("t", good))
expect_error(get_bw(5))
expect_error(get_bw("no_kernel"))
