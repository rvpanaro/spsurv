# Unit tests for print.summary.* (R/print.summary.*.R)

library(spsurv)

context("print.summary")

test_that("print.summary.bpph.mle outputs PH title", {
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle")
  s <- summary(fit)
  expect_output(print(s), "Bernstein Polynomial based Proportional Hazards model")
})

test_that("print.summary.bpph.bayes outputs Bayesian PH title", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1)
  s <- summary(fit)
  expect_output(print(s), "Bayesian Bernstein Polynomial based Proportional Hazards model")
})

test_that("print.summary.bppo.mle outputs PO title", {
  fit <- bppo(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle")
  s <- summary(fit)
  expect_output(print(s), "Bernstein Polynomial based Proportional Odds model")
})

test_that("print.summary.bppo.bayes outputs Bayesian PO title", {
  fit <- bppo(Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1)
  s <- summary(fit)
  expect_output(print(s), "Bayesian Bernstein Polynomial based Proportional Odds model")
})

test_that("print.summary.bpaft.mle outputs AFT title", {
  fit <- bpaft(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle")
  s <- summary(fit)
  expect_output(print(s), "Bernstein Polynomial based Accelerated Failure Time model")
})

test_that("print.summary.bpaft.bayes outputs Bayesian AFT title", {
  fit <- bpaft(Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1)
  s <- summary(fit)
  expect_output(print(s), "Bayesian Bernstein Polynomial based Accelerated Failure Time model")
})

test_that("print.summary.spbp.mle outputs n and coefficients", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle")
  s <- summary(fit)
  expect_output(print(s), "n=")
  expect_output(print(s), "Likelihood ratio test|Wald test")
})

test_that("print.summary.spbp.bayes outputs DIC WAIC LPML", {
  fit <- expect_warning(bpph(Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1))
  s <- summary(fit)
  expect_output(print(s), "DIC=|WAIC=|LPML=")
})
