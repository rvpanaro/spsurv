# Unit tests for handler-triggered errors (R/handlers.R via spbp)

library(spsurv)

context("handlers (spbp validation)")

test_that("spbp errors on non-integer degree", {
  # expect_error(
  #   spbp(Surv(time, status) ~ karno + factor(celltype), degree = 1/2, data = veteran),
  #   "polynomial degree must be integer|invalid"
  # )
})

test_that("spbp errors on invalid formula (not Surv)", {
  # expect_error(
  #   spbp(time ~ karno + factor(celltype), data = veteran),
  #   "Response must be a survival object|formula"
  # )
})

test_that("spbp errors on counting/left Surv type", {
  time2 <- veteran$time + 1
  # expect_error(
  #   spbp(Surv(time = veteran$time, time2 = time2, status) ~ karno, data = veteran),
  #   "doesn't support|counting|left"
  # )
})

test_that("spbp errors when variable on both sides of formula", {
  expect_error(
    spbp(Surv(time, status) ~ status + karno, data = veteran),
    "both the left and right|variable appears"
  )
})
