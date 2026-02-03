# Unit tests for model.matrix.spbp (R/model.matrix.R)

library(spsurv)

context("model.matrix.spbp")

test_that("model.matrix returns features matrix from fit", {
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle")
  mm <- model.matrix(fit)
  expect_equal(mm, fit$features)
  expect_true(is.matrix(mm))
})

test_that("model.matrix has same number of rows as data", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle")
  mm <- model.matrix(fit)
  expect_equal(nrow(mm), nrow(veteran))
})

test_that("model.matrix for null model", {
  fit <- bpph(Surv(time, status) ~ 1, data = veteran, approach = "mle")
  mm <- model.matrix(fit)
  expect_equal(mm, fit$features)
})
