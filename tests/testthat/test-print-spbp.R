# Unit tests for print.spbp (R/print.spbp.R)

library(spsurv)

context("print.spbp")

test_that("print.spbp uses summary printer for MLE fit with covariates", {
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle", init = 0)
  expect_output(print(fit), "Call:")
  expect_output(print(fit), "Bernstein PH model")
  expect_output(print(fit), "Regression coefficients:")
  expect_output(print(fit), "loglik =|AIC =")
})

test_that("print.spbp with bp.param = TRUE", {
  fit <- bpph(Surv(time, status) ~ 1, data = veteran, approach = "mle", init = 0)
  expect_output(print(fit, bp.param = TRUE), "log\\(gamma\\)|gamma")
})

test_that("print.spbp for Bayes fit", {
  fit <- quick_bayes(bpph, Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1, init = 0)
  expect_output(print(fit), "Call:")
  expect_output(print(fit), "Bayesian Bernstein PH model")
  expect_output(print(fit), "Regression coefficients:")
  expect_output(print(fit), "DIC =|WAIC =")
})

test_that("print.spbp for Bayes null model", {
  fit <- quick_bayes(bpph, Surv(time, status) ~ 1, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1, init = 0)
  expect_output(print(fit), "mean\\(bp\\)|mode\\(bp\\)")
})

test_that("print.spbp what = tidy and glance", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  out <- capture.output(print(fit, what = c("tidy", "glance")))
  expect_true(any(grepl("term", out)))
  expect_true(any(grepl("logLik", out)))
  expect_false(any(grepl("waic|dic|lpml", out, ignore.case = TRUE)))
  expect_false(any(grepl("Regression coefficients", out)))
})

test_that("print.spbp glance AIC matches AIC() and logLik nparams", {
  fit <- mle_ph_fit(5L)
  gd <- glance(fit)
  expect_equal(gd$AIC, AIC(fit), tolerance = 1e-10)
  expect_equal(gd$df, attr(logLik(fit), "df"))
})
