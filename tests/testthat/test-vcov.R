# Unit tests for vcov.spbp (R/vcov.R)

library(spsurv)

context("vcov.spbp")

test_that("vcov returns matrix for MLE fit with covariates", {
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle", init = 0)
  v <- vcov(fit)
  expect_true(is.matrix(v))
  expect_equal(dim(v), c(length(coef(fit)), length(coef(fit))))
  expect_equal(rownames(v), names(coef(fit)))
  expect_equal(colnames(v), names(coef(fit)))
})

test_that("vcov with bp.param = TRUE returns full variance matrix", {
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran, approach = "mle", init = 0)
  v <- vcov(fit, bp.param = TRUE)
  expect_true(is.matrix(v))
  q <- length(fit$coefficients)
  m <- length(fit$bp.param)
  expect_equal(dim(v), c(q + m, q + m))
  expect_true(all(is.finite(diag(v)[seq_len(q)])))
})

test_that("vcov with bp.param = TRUE sets gamma variances to NA when ill-conditioned", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
  expect_false(.spbp_gamma_information_diagnostics(fit)$stable)
  expect_silent(v <- vcov(fit, bp.param = TRUE))
  m <- length(fit$bp.param)
  expect_true(all(is.na(diag(v)[(nrow(v) - m + 1L):nrow(v)])))
})

test_that("vcov gives warning and NULL for Bayes fit", {
  fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "bayes", iter = 10, chains = 1, cores = 1)
  expect_warning(v <- vcov(fit), "Not available")
  expect_null(v)
})

.spbp_backtransform_jacobian <- function(object) {
  stan <- object$stanfit
  nms <- names(stan$par)
  eta <- stan$par[startsWith(nms, "beta[")]
  psi <- stan$par[startsWith(nms, "gamma[")]
  means <- object$means
  sdv <- object$sdv
  q <- length(eta)
  m <- length(psi)

  to_orig <- function(e, p) {
    alpha <- sum(e * means / sdv)
    c(beta = e / sdv, gamma = p / exp(alpha))
  }

  eps <- 1e-5
  J <- matrix(0, q + m, q + m)
  for (j in seq_len(q)) {
    e <- eta
    e[j] <- e[j] + eps
    up <- to_orig(e, psi)
    e <- eta
    e[j] <- e[j] - eps
    dn <- to_orig(e, psi)
    J[, j] <- (up - dn) / (2 * eps)
  }
  for (j in seq_len(m)) {
    p <- psi
    p[j] <- p[j] + eps
    up <- to_orig(eta, p)
    p <- psi
    p[j] <- p[j] - eps
    dn <- to_orig(eta, p)
    J[, q + j] <- (up - dn) / (2 * eps)
  }
  J
}

test_that("PH/PO vcov Jacobian matches numerical back-transform", {
  dat <- survival::veteran
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = dat, approach = "mle", init = 0)
  J_num <- .spbp_backtransform_jacobian(fit)

  q <- length(fit$coefficients)
  m <- length(fit$bp.param)
  inv_exp_a <- 1 / exp(fit$alpha)
  J <- diag(q + m)
  diag(J)[seq_len(q)] <- 1 / fit$sdv
  diag(J)[(q + 1):(q + m)] <- inv_exp_a
  for (i in seq_len(m)) {
    J[q + i, seq_len(q)] <- -fit$bp.param[i] * fit$means / fit$sdv
  }

  expect_equal(J, J_num, tolerance = 1e-5)
})

test_that("survfit std.err matches delta method with vcov for PH", {
  dat <- survival::veteran
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = dat, approach = "mle", init = 0)
  t0 <- 50
  sf <- survfit(fit, times = t0)

  q <- length(fit$coefficients)
  m <- length(fit$bp.param)
  X <- t(fit$means)
  G <- sapply(seq_len(m), function(k) pbeta(t0 / fit$tau_b, k, m - k + 1))
  exp_eta <- exp(as.numeric(X %*% fit$coefficients))
  H <- sum(G * fit$bp.param) * exp_eta
  grad <- c(H * as.numeric(X), G * exp_eta)
  V <- vcov(fit, bp.param = TRUE)
  se_manual <- sqrt(as.numeric(t(grad) %*% V %*% grad))

  expect_equal(unname(sf$std.err[1, 1]), unname(se_manual), tolerance = 1e-8)
  expect_lt(se_manual / H, 1)
})

test_that("ill-conditioned gamma block still returns survfit uncertainty with warning", {
  skip_if_not_installed("KMsurv")
  data(larynx, package = "KMsurv")
  larynx$stage <- factor(larynx$stage)
  fit <- bpph(
    Surv(time, delta) ~ age + stage,
    data = larynx,
    approach = "mle",
    init = 0
  )
  expect_false(.spbp_gamma_information_diagnostics(fit)$stable)
  expect_warning(
    sf <- survfit(fit, times = c(10, 50, 100)),
    "lower Bernstein degree"
  )
  expect_true(all(is.finite(sf$std.err)))
  expect_true(all(is.finite(sf$lower)))
  expect_true(all(is.finite(sf$upper)))
  expect_false(any(is.na(sf$surv)))
})

test_that("ill-conditioned gamma block NA only gamma variances in full vcov", {
  skip_if_not_installed("KMsurv")
  data(larynx, package = "KMsurv")
  larynx$stage <- factor(larynx$stage)
  fit <- bpph(
    Surv(time, delta) ~ age + stage,
    data = larynx,
    approach = "mle",
    init = 0
  )
  q <- length(fit$coefficients)
  m <- length(fit$bp.param)
  expect_silent(v_full <- vcov(fit, bp.param = TRUE))
  expect_true(all(is.finite(diag(v_full)[seq_len(q)])))
  expect_true(all(is.na(v_full[(q + 1L):(q + m), (q + 1L):(q + m)])))
  expect_silent(v_beta <- vcov(fit))
  expect_true(all(is.finite(diag(v_beta))))
})

test_that("ill-conditioned gamma warning appears only for survfit, not after fit", {
  skip_if_not_installed("KMsurv")
  data(larynx, package = "KMsurv")
  larynx$stage <- factor(larynx$stage)
  fit <- bpph(
    Surv(time, delta) ~ age + stage,
    data = larynx,
    approach = "mle",
    init = 0
  )
  expect_silent(vcov(fit))
  expect_silent(vcov(fit, bp.param = TRUE))
  expect_silent(summary(fit))
  expect_warning(survfit(fit, times = 50), "lower Bernstein degree")
})
