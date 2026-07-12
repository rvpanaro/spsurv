# Cheap UX wins: rank_models, show_baseline, tidy baseline, ggresiduals, etc.

context("cheap-wins")

test_that("MLE fit does not message about ignored priors by default", {
  fit <- mle_ph_fit(4L)
  expect_silent(fit2 <- bpph(
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "mle",
    degree = 4,
    init = 0,
    verbose = FALSE
  ))
  expect_s3_class(fit2, "spbp")
})

test_that("MLE fit stores residuals by default", {
  fit <- mle_ph_fit(4L)
  expect_true(!is.null(fit$residuals))
  expect_equal(length(fit$residuals), fit$n)
})

test_that("AIC multi-fit table includes model metadata and is sorted", {
  fits <- nested_ph_mle(4L)
  tab <- AIC(null = fits$null, full = fits$full)
  expect_true(all(c("fit", "model", "degree", "aic", "npars") %in% names(tab)))
  expect_equal(tab$aic, sort(tab$aic))
})

test_that("rank_models sorts by AIC", {
  fits <- nested_ph_mle(4L)
  tab <- rank_models(full = fits$full, null = fits$null)
  expect_lt(tab$aic[1L], tab$aic[2L])
})

test_that("summary show_baseline adds baseline table", {
  fit <- mle_ph_fit(4L)
  sm <- summary(fit, show_baseline = TRUE)
  expect_true(!is.null(sm$baseline))
  expect_output(print(sm), "Bernstein baseline coefficients")
})

test_that("tidy baseline component returns gamma rows", {
  fit <- mle_ph_fit(4L)
  td <- tidy(fit, component = "baseline")
  expect_equal(nrow(td), length(fit$bp.param))
  expect_true(all(td$component == "baseline"))
})

test_that("confint works for baseline parameter names", {
  fit <- mle_ph_fit(4L)
  ci <- confint(fit, parm = names(fit$bp.param)[1:2])
  expect_equal(nrow(ci), 2L)
  expect_true(all(ci > 0))
})

test_that("survfit baseline tidy returns time and surv", {
  fit <- mle_ph_fit(4L)
  sf <- survfit(fit, baseline = TRUE, tidy = TRUE)
  expect_true(all(c("time", "surv") %in% names(sf)))
  expect_true(all(sf$surv <= 1))
})

test_that("ggresiduals returns ggplot when ggplot2 available", {
  skip_if_not_installed("ggplot2")
  fit <- mle_ph_fit(4L)
  p <- ggresiduals(fit, type = "martingale")
  expect_s3_class(p, "ggplot")
})
