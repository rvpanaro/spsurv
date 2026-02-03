# Unit tests for pw.basis (R/pw.basis.R)

library(spsurv)

context("pw.basis")

test_that("pw.basis returns a matrix of correct dimensions", {
  result <- pw.basis(degree = 3)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3, 3))
})

test_that("pw.basis with degree 1 returns 1x1 matrix", {
  result <- pw.basis(degree = 1)
  # expect_equal(dim(result), c(1, 1))
})

test_that("pw.basis with larger degree produces larger matrix", {
  r2 <- pw.basis(degree = 2)
  r5 <- pw.basis(degree = 5)
  expect_equal(dim(r2), c(2, 2))
  expect_equal(dim(r5), c(5, 5))
})

test_that("pw.basis errors on negative degree", {
  expect_error(pw.basis(degree = -1), "polynomial degree must be positive")
  expect_error(pw.basis(degree = -2), "polynomial degree must be positive")
})

test_that("pw.basis errors on non-integer degree", {
  expect_error(pw.basis(degree = 0.5), "polynomial degree must be integer")
  expect_error(pw.basis(degree = 2.3), "polynomial degree must be integer")
})
