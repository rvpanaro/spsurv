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
  fit <- quick_bayes(spbp, Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1)
  expect_s3_class(fit, "spbp")
  expect_true(!is.null(fit$posterior$beta) || !is.null(fit$posterior$gamma))
})

test_that("spbp.default uses default degree when missing", {
  fit <- spbp(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  default_degree <- as.integer(ceiling(sqrt(nrow(veteran))))
  expect_equal(fit$degree, default_degree)
  expect_equal(fit$call$degree, default_degree)
  expect_equal(length(fit$bp.param), default_degree)
})

test_that("spbp stores explicit degree in object and call", {
  fit <- spbp(
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "mle",
    degree = 5L,
    init = 0
  )
  expect_equal(fit$degree, 5L)
  expect_equal(fit$call$degree, 5L)
  expect_equal(length(fit$bp.param), 5L)
})

test_that("spbp accepts dist = bernstein(m)", {
  fit <- spbp(
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "mle",
    model = "ph",
    dist = bernstein(4),
    init = 0
  )
  expect_equal(length(fit$bp.param), 4L)
})

test_that("spbp baseline = bernstein(m) equivalent to dist", {
  fit <- spbp(
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "mle",
    model = "ph",
    baseline = bernstein(4),
    init = 0
  )
  expect_equal(length(fit$bp.param), 4L)
})
