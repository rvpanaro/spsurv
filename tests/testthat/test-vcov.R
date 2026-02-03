# Unit tests for vcov.spbp (R/vcov.R)

library(spsurv)

context("vcov.spbp")

test_that("vcov returns matrix for MLE fit with covariates", {
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle")
  v <- vcov(fit)
  expect_true(is.matrix(v))
  expect_equal(dim(v), c(length(coef(fit)), length(coef(fit))))
  expect_equal(rownames(v), names(coef(fit)))
  expect_equal(colnames(v), names(coef(fit)))
})

test_that("vcov with bp.param = TRUE returns full variance matrix", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle")
  v <- vcov(fit, bp.param = TRUE)
  expect_true(is.matrix(v))
  q <- length(fit$coefficients)
  m <- length(fit$bp.param)
  expect_equal(dim(v), c(q + m, q + m))
})

test_that("vcov gives warning and NULL for Bayes fit", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1)
  expect_warning(v <- vcov(fit), "Not available")
  expect_null(v)
})
