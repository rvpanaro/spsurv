# Unit tests for tidy.spbp (R/tidy.spbp.R)

library(spsurv)

context("tidy.spbp")

test_that("tidy.spbp reproduces exp(coef) and exp(confint) workflow for MLE", {
  fit <- bpph(
    Surv(time, status) ~ karno + factor(celltype),
    data = veteran,
    approach = "mle",
    init = 0
  )

  td <- generics::tidy(fit, conf.int = TRUE, exponentiate = TRUE)
  legacy <- cbind(HR = exp(fit$coef), exp(confint(fit)))

  expect_true(all(c(
    "term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high"
  ) %in% names(td)))
  expect_equal(td$term, names(fit$coef))
  expect_equal(unname(td$estimate), unname(legacy[, "HR"]), tolerance = 1e-10)
  expect_equal(unname(td$conf.low), unname(legacy[, 2]), tolerance = 1e-10)
  expect_equal(unname(td$conf.high), unname(legacy[, 3]), tolerance = 1e-10)
})

test_that("tidy.spbp returns expected core columns without intervals", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  td <- generics::tidy(fit)

  expect_true(all(c("term", "estimate", "std.error", "statistic", "p.value") %in% names(td)))
  expect_false(any(c("conf.low", "conf.high") %in% names(td)))
  expect_equal(nrow(td), length(coef(fit)))
})

test_that("tidy.spbp MLE single-coefficient model has non-missing std.error", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  td <- generics::tidy(fit)
  sm <- summary(fit)

  expect_false(any(is.na(td$std.error)))
  expect_equal(td$std.error, sm$coefficients[, "se(coef)"], tolerance = 1e-10)
  expect_equal(td$statistic, sm$coefficients[, "z"], tolerance = 1e-10)
  expect_equal(td$p.value, sm$coefficients[, "Pr(>|z|)"], tolerance = 1e-10)
})

test_that("tidy.spbp Bayes fit returns posterior SD and HPD intervals", {
  fit <- quick_bayes(bpph,
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "bayes",
    iter = 10,
    chains = 1,
    cores = 1,
    init = 0
  )

  td <- generics::tidy(fit, conf.int = TRUE)
  ci <- credint(fit, prob = 0.95, type = "HPD")

  expect_true(all(c("term", "estimate", "std.error", "conf.low", "conf.high") %in% names(td)))
  expect_false(any(c("statistic", "p.value") %in% names(td)))
  expect_false(any(is.na(td$std.error)))
  expect_equal(unname(td$conf.low), unname(ci[, 1]), tolerance = 1e-6)
  expect_equal(unname(td$conf.high), unname(ci[, 2]), tolerance = 1e-6)
})
