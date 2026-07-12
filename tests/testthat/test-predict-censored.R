# Tests for censored-style predict and augment

test_that("predict.spbp survival returns nested censored format", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  pr <- predict(fit, veteran[1:2, ], type = "survival", eval_time = c(100, 200))
  expect_true(".pred" %in% names(pr))
  expect_length(pr$.pred, 2L)
  expect_true(all(c(".eval_time", ".pred_survival") %in% names(pr$.pred[[1L]])))
})

test_that("predict.spbp time returns median survival", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  pr <- predict(fit, veteran[1:3, ], type = "time")
  expect_true(".pred_time" %in% names(pr))
  expect_length(pr$.pred_time, 3L)
})

test_that("predict.spbp linear_pred matches model matrix", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  pr <- predict(fit, veteran[1:3, ], type = "linear_pred")
  X <- model.matrix(fit)[1:3, , drop = FALSE]
  expect_equal(pr$.pred_linear_pred, as.vector(X %*% coef(fit)), tolerance = 1e-10)
})

test_that("predict.spbp curve mode backward compatible", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  pr <- predict(fit, times = c(0, 10, 20))
  expect_true(all(c("time", "surv") %in% names(pr)))
})

test_that("augment.spbp adds residuals", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  aug <- generics::augment(fit, data = veteran[1:5, ])
  expect_true(".residual" %in% names(aug))
  expect_equal(nrow(aug), 5L)
})
