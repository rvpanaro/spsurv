# Unit tests for survfit.spbp and residuals.spbp (R/survfit.R)

library(spsurv)

context("survfit.spbp and residuals.spbp")

test_that("survfit.spbp returns survfit object for MLE", {
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle")
  sf <- survfit(fit)
  expect_s3_class(sf, "survfit")
  expect_true(!is.null(sf$surv))
  expect_true(!is.null(sf$time))
  expect_true(!is.null(sf$cumhaz))
})

test_that("survfit.spbp with newdata", {
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle")
  nd <- veteran[1:2, ]
  sf <- survfit(fit, newdata = nd)
  expect_s3_class(sf, "survfit")
  expect_equal(ncol(sf$surv), 2L)
})

test_that("survfit.spbp for null model", {
  fit <- bpph(Surv(time, status) ~ 1, data = veteran, approach = "mle")
  sf <- survfit(fit)
  expect_s3_class(sf, "survfit")
  expect_true(!is.null(sf$surv))
})

test_that("residuals.spbp returns vector", {
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle")
  r <- residuals(fit)
  expect_true(is.vector(r))
  expect_equal(length(r), nrow(veteran))
  expect_equal(names(r), rownames(fit$features))
})

test_that("residuals.spbp type martingale", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle")
  r <- residuals(fit, type = "martingale")
  expect_true(is.vector(r))
  expect_equal(length(r), nrow(veteran))
})

test_that("residuals.spbp type deviance", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle")
  r <- residuals(fit, type = "deviance")
  expect_true(is.vector(r))
})

test_that("residuals.spbp type cox-snell (coobject-snell in code)", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle")
  r <- residuals(fit, type = "coobject-snell")
  expect_true(is.vector(r))
  expect_true(all(r >= 0))
})
