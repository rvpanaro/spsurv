# Unit tests for summary.spbp (R/summary.spbp.R)

library(spsurv)

context("summary.spbp")

test_that("summary.spbp returns object for MLE fit", {
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle", init = 0)
  s <- summary(fit, interval = 0.95)
  expect_true(is.list(s))
  expect_equal(class(s), "summary.bpph.mle")
  expect_true(!is.null(s$coefficients))
  expect_true(!is.null(s$interval))
  expect_true(!is.null(s$logtest))
  expect_true(!is.null(s$waldtest))
  expect_true(!is.null(s$rsq))
})

test_that("summary.spbp for null model returns fit unchanged", {
  fit <- bpph(Surv(time, status) ~ 1, data = veteran, approach = "mle", init = 0)
  s <- summary(fit)
  expect_equal(s, fit)
})

test_that("summary.spbp for Bayes fit", {
  fit <- expect_warning(bpph(Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1, init = 0))
  s <- summary(fit, interval = 0.95)
  expect_equal(class(s), "summary.bpph.bayes")
  expect_true(!is.null(s$coefficients))
  expect_true(!is.null(s$interval))
  expect_true(!is.null(s$dic) || !is.null(s$waic) || !is.null(s$lpml))
})

test_that("summary.spbp PO model class", {
  fit <- bppo(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  s <- summary(fit)
  expect_equal(class(s), "summary.bppo.mle")
})

test_that("summary.spbp AFT model class", {
  fit <- bpaft(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  s <- summary(fit)
  expect_equal(class(s), "summary.bpaft.mle")
})
