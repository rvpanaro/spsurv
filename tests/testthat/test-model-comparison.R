# Unit tests for bernstein() baseline and model-comparison API (MLE PH/PO/AFT)

library(spsurv)

context("model-comparison")

test_that("bernstein() returns expected structure", {
  b <- bernstein(5)
  expect_equal(b$baseline, "bernstein")
  expect_equal(b$m, 5L)
})

test_that("S3 model-comparison methods are registered for spbp", {
  expect_true(!is.null(getS3method("logLik", "spbp")))
  expect_true(!is.null(getS3method("AIC", "spbp")))
  expect_true(!is.null(getS3method("extractAIC", "spbp")))
  expect_true(!is.null(getS3method("anova", "spbp")))
  expect_true(!is.null(getS3method("estimates", "spbp")))
  expect_true(!is.null(getS3method("se", "spbp")))
})

test_that("logLik.spbp returns logLik object for MLE", {
  fit <- mle_ph_fit(5L)
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), fit$loglik[2], tolerance = 1e-10)
  expect_equal(attr(ll, "df"), length(coef(fit)) + length(fit$bp.param))
  expect_equal(attr(ll, "nobs"), fit$n)
})

test_that("logLik.spbp null model df is length(bp.param)", {
  fits <- nested_ph_mle(5L)
  ll <- logLik(fits$null)
  expect_equal(attr(ll, "df"), length(fits$null$bp.param))
})

test_that("logLik.spbp warns for Bayes fits", {
  fit <- quick_bayes(bpph,
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "bayes",
    iter = 10,
    chains = 1,
    cores = 1,
    init = 0
  )
  expect_warning(res <- logLik(fit))
  expect_null(res)
})

test_that("AIC and BIC use the same nparams as logLik", {
  fit <- mle_ph_fit(5L)
  ll <- logLik(fit)
  k <- attr(ll, "df")
  expect_equal(AIC(fit), -2 * as.numeric(ll) + 2 * k, tolerance = 1e-10)
  gd <- glance(fit)
  expect_equal(gd$df, k)
  expect_equal(gd$AIC, AIC(fit), tolerance = 1e-10)
  expect_equal(gd$BIC, -2 * as.numeric(ll) + log(fit$n) * k, tolerance = 1e-10)
})

test_that("AIC.spbp multi-fit returns comparison data.frame", {
  fits <- nested_ph_mle(5L)
  tab <- AIC(null = fits$null, full = fits$full)
  expect_s3_class(tab, "data.frame")
  expect_true(all(c("fit", "model", "degree", "aic", "npars") %in% names(tab)))
  expect_equal(nrow(tab), 2L)
  expect_equal(tab$npars[tab$fit == "full"], attr(logLik(fits$full), "df"))
})

test_that("extractAIC.spbp returns length-2 vector", {
  fit <- mle_ph_fit(5L)
  ea <- extractAIC(fit)
  expect_length(ea, 2L)
  expect_equal(ea[1L], attr(logLik(fit), "df"))
  expect_equal(ea[2L], AIC(fit), tolerance = 1e-10)
})

test_that("anova.spbp nested PH models", {
  fits <- nested_ph_mle(5L)
  tab <- anova(fits$null, fits$full)
  expect_s3_class(tab, "anova")
  expect_s3_class(tab, "data.frame")
  expect_true(all(c("Terms", "Resid. Df", "-2*LL", "Test", "Df", "Deviance", "Pr(>Chi)") %in% colnames(tab)))
  expect_equal(tab[2L, "Df"], 1L)
  expect_gt(unname(tab[2L, "Deviance"]), 0)
})

test_that("anova.spbp works for PO and AFT nested models", {
  po <- nested_po_mle(5L)
  aft <- nested_aft_mle(5L)
  expect_s3_class(anova(po$null, po$full), "anova")
  expect_s3_class(anova(aft$null, aft$full), "anova")
})

test_that("anova.spbp single fit matches summary logtest", {
  fits <- nested_ph_mle(5L)
  sm <- summary(fits$full)
  tab <- anova(fits$full)
  expect_s3_class(tab, "anova")
  expect_equal(nrow(tab), 2L)
  expect_equal(rownames(tab), c("NULL", "karno"))
  expect_equal(unname(tab["karno", "Deviance"]), unname(sm$logtest["test"]), tolerance = 1e-10)
  expect_equal(unname(tab["karno", "Df"]), unname(sm$logtest["df"]))
})

test_that("anova.spbp sequential table on larynx mirrors survreg layout", {
  skip_if_not_installed("KMsurv")
  data(larynx, package = "KMsurv")
  larynx$stage <- factor(larynx$stage)
  fit0 <- bpph(Surv(time, delta) ~ 1, data = larynx, approach = "mle", degree = 5, init = 0)
  fit1 <- bpph(Surv(time, delta) ~ age, data = larynx, approach = "mle", degree = 5, init = 0)
  fit2 <- bpph(Surv(time, delta) ~ age + stage, data = larynx, approach = "mle", degree = 5, init = 0)

  tab01 <- anova(fit0, fit1)
  tab12 <- anova(fit1, fit2)
  tab2 <- anova(fit2)

  expect_equal(tab01$Terms, c("1", "age"))
  expect_equal(unname(tab01[2L, "Df"]), 1L)
  expect_equal(tab12$Terms, c("age", "age + stage"))
  expect_equal(unname(tab12[2L, "Test"]), "+stage")
  expect_equal(unname(tab12[2L, "Df"]), 3L)

  expect_equal(nrow(tab2), 3L)
  expect_equal(rownames(tab2), c("NULL", "age", "stage"))
  expect_true(is.na(unname(tab2["NULL", "Deviance"])))
  expect_equal(unname(tab2["age", "Df"]), 1L)
  expect_equal(unname(tab2["stage", "Df"]), 3L)
  expect_equal(unname(tab2["age", "Deviance"]), unname(tab01[2L, "Deviance"]), tolerance = 1e-6)
  expect_equal(unname(tab2["stage", "Deviance"]), unname(tab12[2L, "Deviance"]), tolerance = 1e-6)
})

test_that("bpph stores a clean user-facing call", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", degree = 5, init = 0)
  expect_identical(fit$call[[1L]], quote(bpph))
  expect_null(fit$call$dist)
  expect_null(fit$call$baseline)
  expect_equal(fit$call$model, "ph")
})

test_that("anova warns on mismatched Bernstein degree", {
  f0 <- bpph(Surv(time, status) ~ 1, data = veteran, approach = "mle", degree = 4, init = 0)
  f1 <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", degree = 5, init = 0)
  expect_warning(anova(f0, f1))
})

test_that("print.anova works on spbp anova tables", {
  fits <- nested_ph_mle(5L)
  expect_output(print(anova(fits$null, fits$full)), "Analysis of Deviance Table")
})

test_that("estimates and se align for MLE", {
  fit <- mle_ph_fit(5L)
  est <- estimates(fit)
  se_vals <- se(fit)
  expect_equal(length(est), length(se_vals))
  expect_equal(names(se_vals), names(est))
  expect_true(all(names(coef(fit)) %in% names(est)))
  expect_false(any(is.na(se_vals[names(coef(fit))])))
})

test_that("estimates null model length equals bp.param", {
  fits <- nested_ph_mle(5L)
  expect_equal(length(estimates(fits$null)), length(fits$null$bp.param))
})

test_that("spbp resolves dist and baseline to degree", {
  fit_dist <- spbp(
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "mle",
    model = "ph",
    dist = bernstein(4),
    init = 0
  )
  fit_base <- spbp(
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "mle",
    model = "ph",
    baseline = bernstein(4),
    init = 0
  )
  expect_equal(length(fit_dist$bp.param), 4L)
  expect_equal(length(fit_base$bp.param), 4L)
})

test_that("explicit degree overrides dist$m", {
  fit <- bpph(
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "mle",
    degree = 4,
    dist = bernstein(5),
    init = 0
  )
  expect_equal(length(fit$bp.param), 4L)
})

test_that("parametric dist errors with informative message", {
  expect_error(
    bpph(Surv(time, status) ~ karno, data = veteran, dist = "weibull", init = 0),
    "not supported|Parametric"
  )
})

test_that("model-comparison workflow runs end-to-end", {
  fits <- nested_ph_mle(5L)
  expect_s3_class(anova(fits$null, fits$full), "anova")
  expect_s3_class(AIC(fits$null, fits$full), "data.frame")
  td <- tidy(fits$full, conf.int = TRUE)
  expect_true(nrow(td) >= 1L)
  pr <- predict(fits$full, times = c(0, 100, 200))
  expect_true(is.data.frame(pr))
  expect_true(all(c("time", "surv") %in% names(pr)))
})
