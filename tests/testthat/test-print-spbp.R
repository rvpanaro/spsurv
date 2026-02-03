# Unit tests for print.spbp (R/print.spbp.R)

library(spsurv)

context("print.spbp")

test_that("print.spbp runs without error for MLE fit with covariates", {
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle", init = 0)
  expect_output(print(fit), "Call:")
  expect_output(print(fit), "coef")
})

test_that("print.spbp with bp.param = TRUE", {
  fit <- bpph(Surv(time, status) ~ 1, data = veteran, approach = "mle", init = 0)
  expect_output(print(fit, bp.param = TRUE), "log\\(gamma\\)|gamma")
})

test_that("print.spbp for Bayes fit", {
  fit <- expect_warning(bpph(Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1, init = 0))
  expect_output(print(fit), "Call:")
  expect_output(print(fit), "mean\\(coef\\)|Deviance|WAIC")
})

test_that("print.spbp for Bayes null model", {
  fit <- expect_warning(bpph(Surv(time, status) ~ 1, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1, init = 0))
  expect_output(print(fit), "mean\\(bp\\)|mode\\(bp\\)")
})
