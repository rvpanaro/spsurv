# Unit tests for package loading and exports (R/spsurv-package.R)

context("spsurv-package")

test_that("package loads and exports main functions", {
  expect_true("spsurv" %in% loadedNamespaces())
  expect_true(exists("spbp", mode = "function"))
  expect_true(exists("bpph", mode = "function"))
  expect_true(exists("bppo", mode = "function"))
  expect_true(exists("bpaft", mode = "function"))
  expect_true(exists("bp.basis", mode = "function"))
  expect_true(exists("pw.basis", mode = "function"))
})

test_that("veteran data is available", {
  data("cancer", package = "survival")
  expect_true(exists("veteran"))
  expect_true(is.data.frame(veteran))
  expect_true("time" %in% names(veteran))
  expect_true("status" %in% names(veteran))
})
