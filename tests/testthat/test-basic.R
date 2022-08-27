#  File tests/testthat/test-basic.R in package ergm.tapered, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################

expect_summary <- function(s, e, value, coefficients, tolerance=0.001) {
  expect_equal(s, value, tolerance=tolerance, ignore_attr=TRUE)
  expect_equal(coef(e), coefficients, tolerance=tolerance, ignore_attr=TRUE)
}

# a directed nw
data(sampson)
set.seed(42)
c.right <- c(-1.7410903, 0.1456246 )
names(c.right) <- c("edges", "triangle")
test_that("asymmetric, directed", {
  fit <- ergm.tapered(samplike ~ edges + triangles())
  b.sum <- summary(fit)
  expect_summary(b.sum$r, fit, 2, c.right)
})
