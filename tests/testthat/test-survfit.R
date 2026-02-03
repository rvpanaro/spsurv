# Unit tests for survfit.spbp and residuals.spbp (R/survfit.R)

library(spsurv)

context("survfit.spbp and residuals.spbp")

test_that("survfit.spbp returns survfit object for MLE", {
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle", init = 0)
  sf <- survfit(fit)
  expect_s3_class(sf, "survfit")
  expect_true(!is.null(sf$surv))
  expect_true(!is.null(sf$time))
  expect_true(!is.null(sf$cumhaz))
})

test_that("survfit.spbp with newdata", {
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle", init = 0)
  nd <- veteran[1:2, ]
  sf <- survfit(fit, newdata = nd)
  expect_s3_class(sf, "survfit")
  expect_equal(ncol(sf$surv), 2L)
})

test_that("survfit.spbp for null model", {
  fit <- bpph(Surv(time, status) ~ 1, data = veteran, approach = "mle", init = 0)
  sf <- survfit(fit)
  expect_s3_class(sf, "survfit")
  expect_true(!is.null(sf$surv))
})

# Lines 44-54: times argument (must be Surv; else km$time, n.event, n.censor, n.risk updated)
test_that("survfit.spbp errors when times is not a Surv object (lines 44-45)", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  expect_error(survfit(fit, times = 1:10), "times is not a 'Surv' object")
  expect_error(survfit(fit, times = c(1, 2, 3)), "times is not a 'Surv' object")
})

test_that("survfit.spbp with times as Surv object updates time and n.event (lines 46-51)", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  times_surv <- Surv(veteran$time[1:5], veteran$status[1:5])
  sf <- survfit(fit, times = times_surv)
  ord <- order(times_surv[, 1], decreasing = FALSE)
  expect_equal(sf$time, sort(times_surv[, 1], decreasing = FALSE))
  expect_equal(sf$n.event, times_surv[, 2][ord])
  expect_equal(sf$n.censor, 1 - sf$n.event)
  expect_equal(sf$n.risk, length(times_surv[, 1]) - cumsum(sf$n.event) + 1)
})

test_that("residuals.spbp returns vector", {
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle", init = 0)
  r <- residuals(fit)
  expect_true(is.vector(r))
  expect_equal(length(r), nrow(veteran))
  expect_equal(names(r), rownames(fit$features))
})

test_that("residuals.spbp type martingale", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  r <- residuals(fit, type = "martingale")
  expect_true(is.vector(r))
  expect_equal(length(r), nrow(veteran))
})

test_that("residuals.spbp type deviance", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  r <- residuals(fit, type = "deviance")
  expect_true(is.vector(r))
})

test_that("residuals.spbp type cox-snell (coobject-snell in code)", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  r <- residuals(fit, type = "coobject-snell")
  expect_true(is.vector(r))
  expect_true(all(r >= 0))
})
