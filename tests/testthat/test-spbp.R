# Unit tests for spbp and spbp.default (R/spbp.R)

library(spsurv)

context("spbp")

test_that("spbp is generic", {
  expect_true(is.function(spbp))
  expect_true(inherits(spbp, "function"))
})

test_that("spbp.default MLE returns spbp object", {
  fit <- spbp(Surv(time, status) ~ karno + factor(celltype), approach = "mle", data = veteran)
  expect_s3_class(fit, "spbp")
  expect_equal(fit$call$approach, "mle")
  expect_true(!is.null(fit$coefficients) || length(fit$bp.param) > 0L)
})

test_that("spbp.default with model ph, po, aft", {
  fit_ph <- spbp(Surv(time, status) ~ karno, data = veteran, approach = "mle", model = "ph")
  fit_po <- spbp(Surv(time, status) ~ karno, data = veteran, approach = "mle", model = "po")
  fit_aft <- spbp(Surv(time, status) ~ karno, data = veteran, approach = "mle", model = "aft")
  expect_equal(fit_ph$call$model, "ph")
  expect_equal(fit_po$call$model, "po")
  expect_equal(fit_aft$call$model, "aft")
})

test_that("spbp.default Bayes with short run", {
  fit <- spbp(Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1)
  expect_s3_class(fit, "spbp")
  expect_true(!is.null(fit$posterior$beta) || !is.null(fit$posterior$gamma))
})

test_that("spbp.default uses default degree when missing", {
  fit <- spbp(Surv(time, status) ~ karno, data = veteran, approach = "mle")
  expect_true(!is.null(fit$bp.param))
  expect_true(length(fit$bp.param) >= 1L)
})
