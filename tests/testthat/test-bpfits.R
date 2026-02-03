# Unit tests for bpph, bppo, bpaft (R/bpfits.R)

library(spsurv)

context("bpfits: bpph, bppo, bpaft")

test_that("bpph returns spbp object with model ph", {
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle", init = 0)
  expect_s3_class(fit, "spbp")
  expect_equal(fit$call$model, "ph")
  expect_equal(fit$call$approach, "mle")
})

test_that("bpph stores formula and data in call", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, init = 0)
  expect_true(!is.null(fit$call$formula))
  expect_true(!is.null(fit$call$data))
})

test_that("bppo returns spbp object with model po", {
  fit <- bppo(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle", init = 0)
  expect_s3_class(fit, "spbp")
  expect_equal(fit$call$model, "po")
})

test_that("bpaft returns spbp object with model aft", {
  fit <- bpaft(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle", init = 0)
  expect_s3_class(fit, "spbp")
  expect_equal(fit$call$model, "aft")
})

test_that("bpph with approach bayes", {
  fit <- expect_warning(bpph(Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1, init = 0))
  expect_equal(fit$call$approach, "bayes")
  expect_true(!is.null(fit$posterior))
})

test_that("bppo with approach bayes", {
  fit <- expect_warning(bppo(Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1, init = 0))
  expect_equal(fit$call$approach, "bayes")
})

test_that("bpaft with approach bayes", {
  fit <- expect_warning(bpaft(Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1, init = 0))
  expect_equal(fit$call$approach, "bayes")
})

test_that("bpph null model", {
  fit <- bpph(Surv(time, status) ~ 1, data = veteran, approach = "mle")
  expect_s3_class(fit, "spbp")
  expect_true(is.null(fit$coefficients))
})
