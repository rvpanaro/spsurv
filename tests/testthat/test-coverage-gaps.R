# Targeted tests for uncovered paths (covr line coverage).

test_that("bp_survival_reg returns parsnip specs for ph, po, and aft", {
  skip_if_not_installed("parsnip")
  expect_s3_class(spsurv::bp_survival_reg("ph"), "model_spec")
  expect_s3_class(spsurv::bp_survival_reg("po"), "model_spec")
  expect_s3_class(spsurv::bp_survival_reg("aft"), "model_spec")
})

test_that("spsurv_register_parsnip registers proportional_odds engines", {
  skip_if_not_installed("parsnip")
  spsurv_register_parsnip()
  po_eng <- parsnip::show_engines("proportional_odds")$engine
  expect_true("spsurv" %in% po_eng)
  expect_true("spsurv_bayes" %in% po_eng)
  aft_eng <- parsnip::show_engines("survival_reg")$engine
  expect_true("spsurv" %in% aft_eng)
})

test_that("parsnip PO and AFT engines support survival, time, and linear_pred", {
  skip_if_not_installed("parsnip")
  skip_if_not_installed("censored")
  po_spec <- spsurv::proportional_odds() |>
    parsnip::set_engine("spsurv", degree = 4L, scale = FALSE, init = 0)
  po_fit <- parsnip::fit(po_spec, Surv(time, status) ~ karno, data = veteran)
  expect_true(".pred" %in% names(predict(po_fit, veteran[1:2, ], type = "survival", eval_time = 100)))
  expect_true(".pred_time" %in% names(predict(po_fit, veteran[1:2, ], type = "time")))
  expect_true(".pred_linear_pred" %in% names(predict(po_fit, veteran[1:2, ], type = "linear_pred")))

  aft_spec <- parsnip::survival_reg() |>
    parsnip::set_engine("spsurv", degree = 4L, scale = FALSE, init = 0)
  aft_fit <- parsnip::fit(aft_spec, Surv(time, status) ~ karno, data = veteran)
  expect_true(".pred_time" %in% names(predict(aft_fit, veteran[1:2, ], type = "time")))
})

test_that("spsurv parsnip fit warns on case weights", {
  skip_if_not_installed("parsnip")
  expect_warning(
    spsurv::spsurv_fit_proportional_hazards(
      Surv(time, status) ~ karno,
      data = veteran,
      case_weights = rep(1, nrow(veteran)),
      degree = 4L,
      init = 0
    ),
    "Case weights"
  )
})

test_that("spsurv_pred helpers reject non-spbp objects", {
  expect_error(spsurv::spsurv_pred_survival(list(), veteran, 100), "spbp")
  expect_error(spsurv::spsurv_pred_time(list(), veteran), "spbp")
  expect_error(spsurv::spsurv_pred_linear_pred(list(), veteran), "spbp")
})

test_that("spsurv_pred helpers work on wrapped parsnip fits", {
  skip_if_not_installed("parsnip")
  skip_if_not_installed("censored")
  spec <- parsnip::proportional_hazards() |>
    parsnip::set_engine("spsurv", degree = 4L, scale = FALSE, init = 0)
  fit <- parsnip::fit(spec, Surv(time, status) ~ karno, data = veteran)
  pr <- spsurv::spsurv_pred_survival(fit, veteran[1:2, ], eval_time = c(100, 200))
  expect_true(".pred" %in% names(pr))
  pr_t <- spsurv::spsurv_pred_time(fit, veteran[1:2, ])
  expect_true(".pred_time" %in% names(pr_t))
})

test_that("estimates and se cover Bayes gamma and error paths", {
  fit <- quick_bayes(
    bpph,
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "bayes",
    iter = 10L,
    chains = 1L,
    cores = 1L,
    init = 0
  )
  est <- estimates(fit)
  se_vals <- se(fit)
  expect_gt(length(est), length(coef(fit)))
  expect_equal(names(se_vals), names(est))
  expect_error(spsurv:::estimates.spbp(list()), "spbp")
  expect_error(spsurv:::se.spbp(list()), "spbp")

  null_bayes <- quick_bayes(
    bpph,
    Surv(time, status) ~ 1,
    data = veteran,
    approach = "bayes",
    iter = 10L,
    chains = 1L,
    cores = 1L,
    init = 0
  )
  expect_gte(length(estimates(null_bayes)), length(null_bayes$bp.param))
})

test_that("estimates names unnamed gamma for MLE", {
  fit <- bpph(Surv(time, status) ~ 1, data = veteran, approach = "mle", degree = 4L, init = 0)
  est <- estimates(fit)
  expect_true(any(grepl("^gamma", names(est))))
})

test_that("as_draws_df includes gamma and log_lik when requested", {
  skip_if_not_installed("posterior")
  fit <- quick_bayes(
    bpph,
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "bayes",
    iter = 10L,
    chains = 1L,
    cores = 1L,
    init = 0
  )
  dr <- as_draws_df.spbp(fit, variables = c("beta", "gamma", "log_lik"))
  expect_true(any(grepl("^gamma", colnames(dr))))
  expect_true(any(grepl("log_lik", colnames(dr))))
  expect_error(as_draws_df.spbp(
    bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  ), "Bayesian")
})

test_that("gather_draws and spread_draws dispatch on spbp", {
  skip_if_not_installed("posterior")
  skip_if_not_installed("tidybayes")
  fit <- quick_bayes(
    bpph,
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "bayes",
    iter = 10L,
    chains = 1L,
    cores = 1L,
    init = 0
  )
  long <- tidybayes::gather_draws(fit, `beta[karno]`)
  expect_true("beta[karno]" %in% long$.variable)
  wide <- tidybayes::spread_draws(fit, `beta[karno]`)
  expect_true("beta[karno]" %in% names(wide))
})

test_that("spread_surv_draws works for PO and AFT Bayes fits", {
  po_fit <- quick_bayes(
    bppo,
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "bayes",
    iter = 10L,
    chains = 1L,
    cores = 1L,
    init = 0
  )
  po_long <- spread_surv_draws.spbp(po_fit, times = c(50, 100), newdata = veteran[1:2, ])
  expect_true(all(c("time", "surv", ".draw") %in% names(po_long)))

  aft_fit <- quick_bayes(
    bpaft,
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "bayes",
    iter = 10L,
    chains = 1L,
    cores = 1L,
    init = 0
  )
  aft_long <- spread_surv_draws.spbp(aft_fit, times = c(50, 100))
  expect_true(nrow(aft_long) > 0L)
})

test_that("augment supports eval_time and row-name matching", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  aug <- augment(fit, data = veteran[1:5, ], eval_time = c(100, 200))
  expect_true(".pred" %in% names(aug))
  rownames(veteran)[1:5] <- paste0("r", seq_len(5))
  fit2 <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  aug2 <- augment(fit2, data = veteran[1:5, ])
  expect_true(".residual" %in% names(aug2))
})

test_that("predict covers curve conf types and newdata curves", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  pr_loglog <- predict(fit, times = c(0, 50, 100), type = "log-log")
  expect_true(is.data.frame(pr_loglog))
  pr_plain <- predict(fit, veteran[1:2, ], times = seq(0, 200, length.out = 10), type = "plain")
  expect_true(nrow(pr_plain) > 0L)
  expect_error(predict(fit, type = "survival"), "eval_time")
  expect_error(predict(fit, type = "invalid_type"), "type")
})

test_that("predict linear_pred for null MLE model", {
  fit <- bpph(Surv(time, status) ~ 1, data = veteran, approach = "mle", degree = 4L, init = 0)
  pr <- predict(fit, veteran[1:3, ], type = "linear_pred")
  expect_equal(pr$.pred_linear_pred, rep(0, 3L))
})

test_that("ggresiduals supports index axis and coxsnell alias", {
  skip_if_not_installed("ggplot2")
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  p1 <- ggresiduals(fit, type = "martingale", against = "index")
  expect_s3_class(p1, "ggplot")
  p2 <- ggresiduals(fit, type = "coxsnell")
  expect_s3_class(p2, "ggplot")
  expect_error(ggresiduals(list()), "spbp")
})

test_that("confint supports numeric parm and gamma parameters", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  ci_num <- confint(fit, parm = 1L)
  expect_true(is.matrix(ci_num))
  gname <- names(fit$bp.param)[1L]
  ci_g <- confint(fit, parm = gname, level = 0.9)
  expect_true(is.matrix(ci_g))
  expect_error(confint(fit, parm = "not_a_param"), "No matching")
})

test_that("credint supports null Bayes gamma intervals", {
  fit <- quick_bayes(
    bpph,
    Surv(time, status) ~ 1,
    data = veteran,
    approach = "bayes",
    iter = 10L,
    chains = 1L,
    cores = 1L,
    init = 0
  )
  ci_hpd <- credint(fit, prob = 0.95, type = "HPD")
  expect_true(is.matrix(ci_hpd))
  ci_eq <- credint(fit, prob = 0.95, type = "Equal-Tailed")
  expect_true(is.matrix(ci_eq))
})

test_that("glance handles degenerate spbp objects", {
  empty <- structure(
    list(
      n = 10L,
      call = list(approach = "mle", model = "ph"),
      coefficients = NULL,
      bp.param = NULL
    ),
    class = "spbp"
  )
  gd <- glance(empty)
  expect_equal(nrow(gd), 1L)
  expect_equal(gd$n, 10L)
})

test_that("tidy supports baseline component and exponentiate without confint", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  base <- tidy(fit, component = "baseline")
  expect_true(all(c("term", "component", "estimate") %in% names(base)))
  td <- tidy(fit, exponentiate = TRUE, conf.int = FALSE)
  expect_true(all(td$estimate > 0))
  td_coef <- tidy(fit, component = "coef", conf.int = TRUE)
  expect_gte(nrow(td_coef), 1L)
})

test_that("tidy null coefficient model returns empty coef frame", {
  fit <- bpph(Surv(time, status) ~ 1, data = veteran, approach = "mle", degree = 4L, init = 0)
  td <- tidy(fit, component = "coef")
  expect_equal(nrow(td), 0L)
})

test_that("coef Bayes null posterior returns NULL", {
  fit <- quick_bayes(
    bpph,
    Surv(time, status) ~ 1,
    data = veteran,
    approach = "bayes",
    iter = 10L,
    chains = 1L,
    cores = 1L,
    init = 0
  )
  fit$posterior$beta <- NULL
  expect_null(coef(fit))
})

test_that("residuals cover PO, AFT, and Bayes paths", {
  po_fit <- bppo(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  expect_length(residuals(po_fit, type = "deviance"), nrow(veteran))
  aft_fit <- bpaft(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  expect_length(residuals(aft_fit, type = "martingale"), nrow(veteran))
  bayes_fit <- quick_bayes(
    bppo,
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "bayes",
    iter = 10L,
    chains = 1L,
    cores = 1L,
    init = 0
  )
  expect_length(residuals(bayes_fit, type = "coxsnell"), nrow(veteran))
  null_aft <- bpaft(Surv(time, status) ~ 1, data = veteran, approach = "mle", degree = 4L, init = 0)
  expect_length(residuals(null_aft), nrow(veteran))
})

test_that("survfit as.data.frame and tidy newdata merge", {
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle", init = 0)
  sf <- survfit(fit, times = seq(0, 200, length.out = 20))
  df <- as.data.frame(sf)
  expect_true(all(c("id", "time", "surv", "lower", "upper") %in% names(df)))
  nd <- veteran[1:2, ]
  sf2 <- survfit(fit, newdata = nd, times = seq(0, 100, by = 10), tidy = TRUE)
  expect_true(any(names(sf2) %in% names(nd)))
})

test_that("internal utils helpers are exercised", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  expect_true(is.data.frame(spsurv:::.spbp_training_data(fit)))
  expect_gt(length(spsurv:::.spbp_default_survfit_times(fit, length.out = 50L)), 1L)
  bands <- spsurv:::.spbp_monotone_surv_bands(
    lower = matrix(c(0.9, 0.8, 0.7), ncol = 1L),
    upper = matrix(c(1.0, 0.95, 0.85), ncol = 1L)
  )
  expect_true(all(bands$upper >= bands$lower))
  expect_type(spsurv:::.spbp_default_cores(), "integer")
  bayes_fit <- quick_bayes(
    bpph,
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "bayes",
    iter = 10L,
    chains = 1L,
    cores = 1L,
    init = 0
  )
  idx <- spsurv:::.spbp_posterior_draw_indices(
    bayes_fit$stanfit,
    nrow(bayes_fit$posterior$beta)
  )
  expect_true(all(c(".chain", ".iteration", ".draw") %in% names(idx)))
  iv <- spsurv:::.spbp_posterior_interval(rnorm(100), interval = 0.9, interval.type = "quantile")
  expect_length(iv, 2L)
})

test_that("print.summary supports baseline table and interval controls", {
  fit <- bpph(Surv(time, status) ~ 1, data = veteran, approach = "mle", degree = 4L, init = 0)
  s <- summary(fit, show_intervals = TRUE, show_baseline = TRUE, compact = FALSE)
  expect_output(print(s), "gamma|log\\(gamma\\)")
  fit2 <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  s2 <- summary(fit2, show_intervals = FALSE)
  expect_output(print(s2), "Regression coefficients")
})

test_that("print.spbp what = glance for Bayes fit", {
  fit <- quick_bayes(
    bpph,
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "bayes",
    iter = 10L,
    chains = 1L,
    cores = 1L,
    init = 0
  )
  out <- capture.output(print(fit, what = "glance"))
  expect_true(any(grepl("waic|dic|lpml", out, ignore.case = TRUE)))
})

test_that("spsurv_register_tidybayes_s3 registers posterior hooks", {
  skip_if_not_installed("posterior")
  spsurv_register_tidybayes_s3()
  expect_true(is.function(
    getS3method("as_draws_df", "spbp", envir = asNamespace("posterior"), optional = TRUE)
  ))
})

test_that("survfit baseline curve supports se.fit and Surv times", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0, degree = 4L)
  sf <- survfit(fit, baseline = TRUE, se.fit = TRUE, times = c(0, 50, 100))
  expect_s3_class(sf, "survfit")
  expect_true(all(c("time", "surv", "std.err", "lower", "upper") %in% names(sf)))
  sf2 <- survfit(
    fit,
    baseline = TRUE,
    times = Surv(veteran$time[1:5], veteran$status[1:5])
  )
  expect_length(sf2$time, 5L)
})

test_that("spbp resolves cores when cores is invalid", {
  fit <- bpph(
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "mle",
    cores = NA,
    degree = 4L,
    init = 0
  )
  expect_s3_class(fit, "spbp")
})

test_that("spbp degree resolution warns on both dist and baseline", {
  expect_warning(
    fit <- bpph(
      Surv(time, status) ~ karno,
      data = veteran,
      approach = "mle",
      dist = bernstein(5L),
      baseline = bernstein(6L),
      init = 0
    ),
    "Both 'dist' and 'baseline'"
  )
  expect_equal(length(fit$bp.param), 5L)
})

test_that("explicit degree overrides dist and baseline without warning", {
  fit <- bpph(
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "mle",
    degree = 4L,
    dist = bernstein(5L),
    baseline = bernstein(6L),
    init = 0
  )
  expect_equal(length(fit$bp.param), 4L)
})

test_that("predict survival uses training data when newdata is NULL", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0, degree = 4L)
  pr <- predict(fit, newdata = NULL, type = "survival", eval_time = c(50, 100))
  expect_length(pr$.pred[[1L]]$.pred_survival, 2L)
})

test_that("internal predict helpers handle edge cases", {
  expect_equal(
    spsurv:::.spbp_interp_surv(c(1, 2), c(0.9, 0.8), c(1.5, 3)),
    stats::approx(c(1, 2), c(0.9, 0.8), xout = c(1.5, 3), rule = 2, ties = mean)$y,
    tolerance = 1e-10
  )
  expect_true(all(is.na(spsurv:::.spbp_interp_surv(c(NA, 1), c(0.5, NA), 1:3))))
  expect_equal(spsurv:::.spbp_interp_surv(5, 0.7, 1:3), rep(0.7, 3))
  expect_error(spsurv:::.spbp_resolve_predict_mode("bad"), "type")
  expect_equal(spsurv:::.spbp_resolve_predict_mode("curve")$mode, "curve")
})

test_that("confint returns gamma-only intervals on log scale", {
  fit <- bpph(Surv(time, status) ~ 1, data = veteran, approach = "mle", degree = 4L, init = 0)
  g0 <- names(fit$bp.param)[1L]
  ci <- confint(fit, parm = g0)
  expect_true(is.matrix(ci))
  expect_true(all(ci > 0))
})

test_that("credint and confint warn on mismatched approach for gamma", {
  fit <- quick_bayes(
    bpph,
    Surv(time, status) ~ 1,
    data = veteran,
    approach = "bayes",
    iter = 10L,
    chains = 1L,
    cores = 1L,
    init = 0
  )
  expect_warning(confint(fit, parm = names(fit$bp.param)[1L]), "credint")
})

test_that("Bayes AFT residuals path is covered", {
  aft_fit <- quick_bayes(
    bpaft,
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "bayes",
    iter = 10L,
    chains = 1L,
    cores = 1L,
    init = 0
  )
  expect_length(residuals(aft_fit, type = "deviance"), nrow(veteran))
})

test_that("as_draws_df synthesizes draw indices when missing", {
  skip_if_not_installed("posterior")
  fit <- quick_bayes(
    bpph,
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "bayes",
    iter = 10L,
    chains = 1L,
    cores = 1L,
    init = 0
  )
  fit$posterior$draw_indices <- NULL
  dr <- as_draws_df.spbp(fit, variables = "beta")
  expect_true(all(c(".chain", ".iteration", ".draw") %in% names(dr)))
})

test_that("gamma information diagnostics handle missing Hessian blocks", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0, degree = 4L)
  diag <- spsurv:::.spbp_gamma_information_diagnostics(fit)
  expect_type(diag$stable, "logical")
  fit$hessian <- NULL
  hx <- spsurv:::.spbp_hessian_neg_beta_gamma(fit)
  expect_null(hx$H_neg)
})

test_that("print.spbp interval display paths", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0, degree = 4L)
  expect_output(print(fit, interval = 0.95), "2.5 %|2.5%")
  bayes_fit <- quick_bayes(
    bpph,
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "bayes",
    iter = 10L,
    chains = 1L,
    cores = 1L,
    init = 0
  )
  expect_output(print(bayes_fit, interval = 0.95), "2.5 %|2.5%")
})
