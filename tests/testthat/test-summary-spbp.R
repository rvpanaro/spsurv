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
  expect_true(!is.null(s$coef_interval))
  expect_true(!is.null(s$logtest))
  expect_true(!is.null(s$waldtest))
  expect_true(!is.null(s$rsq))
  expect_equal(s$nparams, length(coef(fit)) + length(fit$bp.param))
})

test_that("summary.spbp for null model returns fit unchanged", {
  fit <- bpph(Surv(time, status) ~ 1, data = veteran, approach = "mle", init = 0)
  s <- summary(fit)
  expect_equal(s, fit)
})

test_that("summary.spbp for Bayes fit", {
  fit <- quick_bayes(bpph, Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1, init = 0)
  s <- summary(fit, interval = 0.95)
  expect_equal(class(s), "summary.bpph.bayes")
  expect_true(!is.null(s$coefficients))
  expect_true(!is.null(s$interval))
  expect_true(!is.null(s$coef_interval))
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

test_that("summary.spbp coef_interval matches Wald bounds for MLE", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  s <- summary(fit, interval = 0.95)
  z <- qnorm(0.975)
  se <- s$coefficients[, "se(coef)"]
  expect_equal(
    s$coef_interval[, "lower .95"],
    s$coef_interval[, "coef"] - z * se,
    tolerance = 1e-10
  )
  expect_equal(
    s$coef_interval[, "upper .95"],
    s$coef_interval[, "coef"] + z * se,
    tolerance = 1e-10
  )
})

test_that("print.summary hides coef intervals when show_intervals = FALSE", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  s <- summary(fit, show_intervals = FALSE)
  out <- capture.output(print(s))
  expect_false(any(grepl("2\\.5%", out)))
})

test_that("summary logtest consistent with logLik and nested anova", {
  fits <- nested_ph_mle(5L)
  sm <- summary(fits$full)
  expect_equal(
    unname(sm$logtest["test"]),
    -2 * (fits$full$loglik[1] - fits$full$loglik[2]),
    tolerance = 1e-6
  )
  expect_equal(unname(sm$logtest["df"]), length(coef(fits$full)))
})
