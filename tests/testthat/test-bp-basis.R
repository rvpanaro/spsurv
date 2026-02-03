# Unit tests for bp.basis (R/bp.basis.R)

library(spsurv)

context("bp.basis")

test_that("bp.basis returns correct structure", {
  time <- c(1, 2, 3, 4, 5)
  degree <- 3
  result <- bp.basis(time, degree)

  expect_type(result, "list")
  expect_named(result, c("g", "G", "degree", "tau"))
  expect_equal(result$degree, degree)
  expect_equal(result$tau, max(time))

  expect_true(is.matrix(result$g))
  expect_true(is.matrix(result$G))
  expect_equal(nrow(result$g), length(time))
  expect_equal(ncol(result$g), degree)
  expect_equal(dim(result$g), dim(result$G))
})

test_that("bp.basis uses default tau when omitted", {
  time <- c(0.5, 1, 2, 3)
  degree <- 2
  result <- bp.basis(time, degree)

  expect_equal(result$tau, 3)
  expect_equal(result$tau, max(time))
})

test_that("bp.basis accepts explicit tau", {
  time <- c(1, 2, 3)
  degree <- 2
  tau <- 10
  result <- bp.basis(time, degree, tau = tau)

  expect_equal(result$tau, tau)
})

# test_that("bp.basis works with single time value", {
#   result <- bp.basis(time = 1, degree = 2)
#
#   expect_equal(nrow(result$g), 1)
#   expect_equal(ncol(result$g), 2)
#   expect_equal(nrow(result$G), 1)
#   expect_equal(ncol(result$G), 2)
# })

test_that("bp.basis works with degree 1", {
  time <- c(1, 2, 3)
  result <- bp.basis(time, degree = 1)

  expect_equal(dim(result$g), c(3, 1))
  expect_equal(dim(result$G), c(3, 1))
  expect_equal(result$degree, 1)
})

test_that("bp.basis G (CDF) is in [0, 1]", {
  time <- seq(0.1, 5, length.out = 10)
  result <- bp.basis(time, degree = 4, tau = 6)

  expect_true(all(result$G >= 0 & result$G <= 1))
})

test_that("bp.basis g (density) is non-negative", {
  time <- c(0.5, 1, 2, 3, 4)
  result <- bp.basis(time, degree = 3, tau = 5)

  expect_true(all(result$g >= 0))
})

# test_that("bp.basis at time = 0 gives G = 0 for first basis", {
#   result <- bp.basis(time = 0.001, degree = 2, tau = 1)
#   # pbeta(0, k, degree-k+1) = 0 for any k
#   expect_true(all(result$G[1, ] >= 0))
# })

# test_that("bp.basis at time = tau gives G = 1", {
#   tau <- 5
#   result <- bp.basis(time = tau, degree = 3, tau = tau)
#   # pbeta(1, k, m-k+1) = 1
#   expect_equal(result$G[1, ], c(1, 1, 1))
# })

test_that("bp.basis errors on negative time", {
  expect_error(bp.basis(time = -1, degree = 1), "time must be a positive vector")
  expect_error(bp.basis(time = c(1, 2, -3), degree = 2), "time must be a positive vector")
})

test_that("bp.basis errors on zero time when it implies non-positive", {
  # time >= 0 is required: sum(time >= 0) == n. So 0 is allowed.
  result <- bp.basis(time = 0, degree = 1, tau = 1)
  expect_equal(result$tau, 1)
})

test_that("bp.basis errors on negative degree", {
  expect_error(bp.basis(time = 1, degree = -1), "polynomial degree must be positive")
  expect_error(bp.basis(time = 1, degree = -2), "polynomial degree must be positive")
})

test_that("bp.basis errors on non-integer degree", {
  expect_error(bp.basis(time = 1, degree = 0.5), "polynomial degree must be integer")
  expect_error(bp.basis(time = 1, degree = 2.3), "polynomial degree must be integer")
})

test_that("bp.basis errors when tau is less than max(time)", {
  expect_error(bp.basis(time = c(1, 2, 5), degree = 1, tau = 3), "tau must be greater than the last time")
  expect_error(bp.basis(time = 10, degree = 1, tau = 5), "tau must be greater than the last time")
})

test_that("bp.basis accepts tau equal to max(time)", {
  time <- c(1, 2, 3)
  result <- bp.basis(time, degree = 2, tau = 3)
  expect_equal(result$tau, 3)
})

test_that("bp.basis g and G have consistent dimensions with input", {
  time <- c(0.5, 1, 1.5, 2, 2.5)
  degree <- 5
  result <- bp.basis(time, degree)

  expect_equal(nrow(result$g), 5)
  expect_equal(ncol(result$g), 5)
  expect_equal(nrow(result$G), 5)
  expect_equal(ncol(result$G), 5)
})

test_that("bp.basis with larger degree produces more columns", {
  time <- c(1, 2, 3)
  r2 <- bp.basis(time, degree = 2)
  r4 <- bp.basis(time, degree = 4)

  expect_equal(ncol(r2$g), 2)
  expect_equal(ncol(r4$g), 4)
  expect_equal(ncol(r2$G), 2)
  expect_equal(ncol(r4$G), 4)
})
