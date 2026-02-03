# Unit tests for coef.spbp (R/coef.R)

library(spsurv)

context("coef.spbp")

test_that("coef returns coefficients for MLE fit", {
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle")
  b <- coef(fit)
  # expect_true(is.vector(b))
  expect_equal(names(b), names(fit$coefficients))
  expect_equal(b, fit$coefficients)
})

test_that("coef returns coefficients for MLE null model", {
  fit <- bpph(Surv(time, status) ~ 1, data = veteran, approach = "mle")
  b <- coef(fit)
  expect_true(is.null(b) || length(b) == 0L)
})

test_that("coef with summary mean for Bayes fit (line 19-20)", {
  fit <- expect_warning(bpph(Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1))
  b <- coef(fit, summary = "mean")
  expect_true(is.vector(b))
  expect_equal(length(b), ncol(fit$posterior$beta))
  expect_equal(b, apply(fit$posterior$beta, 2, mean))
})

test_that("coef with summary median for Bayes fit (line 21-22)", {
  fit <- expect_warning(bpph(Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1))
  b <- coef(fit, summary = "median")
  expect_true(is.vector(b))
  expect_equal(length(b), ncol(fit$posterior$beta))
  expect_true(is.numeric(b))
})

test_that("coef with summary mode for Bayes fit (line 23-24)", {
  fit <- expect_warning(bpph(Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1))
  b <- coef(fit, summary = "mode")
  expect_true(is.vector(b))
  expect_equal(length(b), ncol(fit$posterior$beta))
  expect_true(is.numeric(b))
})

test_that("coef default summary is mean for Bayes fit (line 16-17)", {
  fit <- expect_warning(bpph(Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1))
  expect_equal(coef(fit), coef(fit, summary = "mean"))
})
