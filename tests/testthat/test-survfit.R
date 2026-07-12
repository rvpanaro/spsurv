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
  expect_gt(length(sf$time), 50L)
})

test_that("bayesian survfit default grid yields monotone credible bands", {
  fit <- expect_warning(
    bpph(
      Surv(time, status) ~ karno,
      data = veteran,
      approach = "bayes",
      iter = 40,
      chains = 1,
      cores = 1,
      init = 0
    )
  )
  sf <- survfit(fit, type = "plain", interval.type = "hpd")
  expect_gt(length(sf$time), 50L)
  lower <- sf$lower[, 1]
  upper <- sf$upper[, 1]
  expect_true(all(diff(lower) <= sqrt(.Machine$double.eps) * pmax(1, abs(lower[-1]))))
  expect_true(all(diff(upper) <= sqrt(.Machine$double.eps) * pmax(1, abs(upper[-1]))))
  expect_true(all(upper >= lower))
})

test_that("survfit.spbp with newdata", {
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle", init = 0)
  nd <- veteran[1:2, ]
  sf <- survfit(fit, newdata = nd)
  expect_s3_class(sf, "survfit")
  expect_equal(ncol(sf$surv), 2L)
})

test_that("survfit.spbp resolves training data via call stack when data slot is absent", {
  local({
    dat <- survival::veteran
    fit <- bpph(Surv(time, status) ~ karno, data = dat, approach = "mle", init = 0)
    fit$data <- NULL
    sf <- survfit(fit)
    expect_s3_class(sf, "survfit")
  })
})

test_that("survfit.spbp for null model", {
  fit <- bpph(Surv(time, status) ~ 1, data = veteran, approach = "mle", init = 0)
  sf <- survfit(fit)
  expect_s3_class(sf, "survfit")
  expect_true(!is.null(sf$surv))
})

test_that("survfit.spbp errors when times is neither Surv nor numeric", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  expect_error(survfit(fit, times = letters[1:3]), "'times' must be numeric or a Surv object")
})

test_that("survfit.spbp accepts numeric times for smooth evaluation grid", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  g <- seq(0, 200, by = 5)
  sf <- survfit(fit, times = g)
  expect_equal(sf$time, g)
  expect_equal(nrow(sf$surv), length(g))
})

test_that("as.data.frame.survfitbp stacks curves for ggplot", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  nd <- veteran[1:2, ]
  sf <- survfit(fit, newdata = nd, times = c(0, 10, 20, 50))
  df <- as.data.frame(sf)
  expect_equal(nrow(df), nrow(sf$surv) * ncol(sf$surv))
  expect_true(all(c("id", "time", "surv", "lower", "upper") %in% names(df)))
})

test_that("predict.spbp returns tidy survival data frame", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  pr <- predict(fit, times = seq(0, 100, length.out = 25))
  expect_true(is.data.frame(pr))
  expect_equal(nrow(pr), 25L)
  expect_true(all(pr$time >= 0))
})

test_that("bpaft predict returns surv = 1 at time = 0", {
  fit <- bpaft(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  pr <- predict(fit, times = c(0, 50, 100))
  expect_equal(pr$surv[pr$time == 0], 1)
  expect_false(any(is.nan(pr$surv)))
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

test_that("residuals.spbp type cox-snell", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  r <- residuals(fit, type = "cox-snell")
  expect_true(is.vector(r))
  expect_true(all(r >= 0))
})

test_that("survfit.spbp tidy = TRUE returns data.frame", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  pr <- survfit(fit, times = c(0, 10, 20, 50), tidy = TRUE)
  expect_true(is.data.frame(pr))
  expect_true(all(c("time", "surv") %in% names(pr)))
})

test_that("survfit tidy with newdata includes covariate columns", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  nd <- veteran[1:2, , drop = FALSE]
  g <- c(0, 10, 20)
  pr <- survfit(fit, newdata = nd, times = g, tidy = TRUE)
  expect_true("karno" %in% names(pr))
  expect_equal(nrow(pr), length(g) * nrow(nd))
})

test_that("predict matches survfit tidy for same times", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  g <- c(0, 10, 20, 50)
  pr <- predict(fit, times = g)
  sf <- survfit(fit, times = g, tidy = TRUE)
  expect_equal(pr$surv, sf$surv, tolerance = 1e-10)
})

test_that("residuals coxsnell alias equals cox-snell", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  expect_equal(
    residuals(fit, type = "coxsnell"),
    residuals(fit, type = "cox-snell")
  )
})
