# Unit tests for glance.spbp (R/glance.spbp.R)

library(spsurv)

context("glance.spbp")

test_that("glance.spbp returns one-row MLE summary aligned with summary()", {
  fit <- bpph(
    Surv(time, status) ~ karno + factor(celltype),
    data = veteran,
    approach = "mle",
    init = 0
  )
  sm <- summary(fit)
  gd <- generics::glance(fit)

  expect_equal(nrow(gd), 1L)
  expect_true(all(c(
    "n", "nevent", "logLik", "approach", "model", "df",
    "statistic", "p.value", "rsq", "max.rsq", "AIC", "BIC"
  ) %in% names(gd)))
  expect_equal(gd$n, fit$n)
  expect_equal(gd$nevent, fit$nevent)
  expect_equal(gd$logLik, as.numeric(logLik(fit)), tolerance = 1e-10)
  expect_equal(gd$approach, "mle")
  expect_equal(gd$model, "ph")
  expect_equal(gd$df, attr(logLik(fit), "df"))
  expect_equal(gd$AIC, AIC(fit), tolerance = 1e-10)
  expect_equal(gd$BIC, -2 * as.numeric(logLik(fit)) + log(fit$n) * attr(logLik(fit), "df"),
    tolerance = 1e-10
  )
  expect_equal(gd$statistic, unname(sm$logtest["test"]), tolerance = 1e-10)
  expect_equal(gd$p.value, unname(sm$logtest["pvalue"]), tolerance = 1e-10)
  expect_equal(gd$rsq, unname(sm$rsq["rsq"]), tolerance = 1e-10)
  expect_equal(gd$max.rsq, unname(sm$rsq["maxrsq"]), tolerance = 1e-10)
  expect_false(any(c("waic", "dic", "lpml") %in% names(gd)))
})

test_that("glance.spbp returns one-row Bayes summary aligned with summary()", {
  fit <- quick_bayes(bpph,
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "bayes",
    iter = 10,
    chains = 1,
    cores = 1,
    init = 0
  )
  sm <- summary(fit)
  gd <- generics::glance(fit)

  expect_equal(nrow(gd), 1L)
  expect_equal(gd$approach, "bayes")
  expect_equal(gd$logLik, sum(fit$loglik), tolerance = 1e-10)
  expect_false(any(c(
    "statistic", "p.value", "rsq", "max.rsq", "AIC", "BIC"
  ) %in% names(gd)))
  expect_equal(gd$waic, sm$waic, tolerance = 1e-6)
  expect_equal(gd$dic, sm$dic, tolerance = 1e-6)
  expect_equal(gd$lpml, sm$lpml, tolerance = 1e-6)
})

test_that("glance.spbp null MLE model reports logLik and AIC from bp.param only", {
  fit <- bpph(Surv(time, status) ~ 1, data = veteran, approach = "mle", degree = 5, init = 0)
  gd <- generics::glance(fit)
  ll <- logLik(fit)
  expect_equal(gd$logLik, as.numeric(ll), tolerance = 1e-10)
  expect_equal(gd$df, attr(ll, "df"))
  expect_equal(gd$df, length(fit$bp.param))
  expect_equal(gd$AIC, AIC(fit), tolerance = 1e-10)
  expect_false(any(c("statistic", "p.value", "rsq") %in% names(gd)))
})
