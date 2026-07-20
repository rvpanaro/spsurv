# Unit tests for print.summary.* (R/print.summary.*.R)

library(spsurv)

context("print.summary")

test_that("print.summary.bpph.mle outputs standard PH summary", {
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle", init = 0)
  s <- summary(fit)
  expect_output(print(s), "Call:")
  expect_output(print(s), "Bernstein PH model")
  expect_output(print(s), "Regression coefficients:")
  expect_output(print(s), "2.5%")
  expect_output(print(s), "Exponentiated coefficients:")
  expect_output(print(s), "loglik =|AIC =")
})

test_that("print.summary.bpph.bayes outputs standard Bayesian PH summary", {
  fit <- quick_bayes(bpph, Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1, init = 0)
  s <- summary(fit)
  expect_output(print(s), "Call:")
  expect_output(print(s), "Bayesian Bernstein PH model")
  expect_output(print(s), "Regression coefficients:")
  expect_output(print(s), "2.5%")
  expect_output(print(s), "DIC =|WAIC =")
})

test_that("print.summary.bppo.mle outputs PO title", {
  fit <- bppo(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle", init = 0)
  s <- summary(fit)
  expect_output(print(s), "Bernstein PO model")
})

test_that("print.summary.bppo.bayes outputs Bayesian PO title", {
  fit <- quick_bayes(bppo, Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1, init = 0)
  s <- summary(fit)
  expect_output(print(s), "Bayesian Bernstein PO model")
})

test_that("print.summary.bpaft.mle outputs AFT title", {
  fit <- bpaft(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle", init = 0)
  s <- summary(fit)
  expect_output(print(s), "Bernstein AFT model")
})

test_that("print.summary.bpaft.bayes outputs Bayesian AFT title", {
  fit <- quick_bayes(bpaft, Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1, init = 0)
  s <- summary(fit)
  expect_output(print(s), "Bayesian Bernstein AFT model")
})

test_that("print.summary.spbp.mle compact omits global tests", {
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle", init = 0)
  s <- summary(fit, compact = TRUE)
  out <- capture.output(print(s))
  expect_match(paste(out, collapse = "\n"), "loglik =")
  expect_false(any(grepl("Likelihood ratio test", out)))
})

test_that("print.summary MLE AIC matches AIC()", {
  fit <- mle_ph_fit(5L)
  s <- summary(fit)
  out <- paste(capture.output(print(s)), collapse = "\n")
  aic_val <- AIC(fit)
  aic_line <- regmatches(out, regexpr("AIC = [0-9.]+", out))
  printed_aic <- as.numeric(sub("AIC = ", "", aic_line))
  expect_equal(printed_aic, aic_val, tolerance = 0.5)
  expect_equal(s$nparams, attr(logLik(fit), "df"))
})

test_that("print.summary.spbp.mle non-compact includes global tests", {
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle", init = 0)
  s <- summary(fit, compact = FALSE)
  expect_output(print(s), "Likelihood ratio test|Wald test")
})

test_that("print.summary.spbp.bayes outputs DIC WAIC LPML when not compact", {
  fit <- quick_bayes(bpph, Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1, init = 0)
  s <- summary(fit, compact = FALSE)
  expect_output(print(s), "DIC =|WAIC =|LPML =")
})
