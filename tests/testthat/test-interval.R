# Unit tests for confint.spbp and credint.spbp (R/interval.R)

library(spsurv)

context("confint.spbp and credint.spbp")

test_that("confint returns matrix for MLE fit with covariates", {
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle")
  ci <- confint(fit, level = 0.95)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), length(coef(fit)))
  expect_equal(ncol(ci), 2L)
  expect_true(all(ci[, 1] <= ci[, 2]))
})

test_that("confint with custom level", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle")
  ci <- confint(fit, level = 0.9)
  expect_true(is.matrix(ci))
})

test_that("credint returns matrix for Bayes fit", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1)
  ci <- credint(fit, prob = 0.95, type = "Equal-Tailed")
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), ncol(fit$posterior$beta))
  expect_equal(ncol(ci), 2L)
})

test_that("credint with type HPD", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1)
  ci <- credint(fit, prob = 0.95, type = "HPD")
  expect_true(is.matrix(ci))
})

test_that("confint on Bayes fit warns and calls credint", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1)
  expect_warning(ci <- confint(fit), "not applicable")
  expect_true(is.matrix(ci))
})

test_that("credint on MLE fit warns and calls confint", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle")
  expect_warning(ci <- credint(fit), "not applicable")
  expect_true(is.matrix(ci))
})

test_that("confint for null MLE model (bp.param)", {
  fit <- bpph(Surv(time, status) ~ 1, data = veteran, approach = "mle")
  ci <- confint(fit, level = 0.95)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), length(fit$bp.param))
})
