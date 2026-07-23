# Finite-difference checks for survfit variance gradients

library(spsurv)

context("survfit variance gradients")

.fd_grad <- function(f, x, eps = 1e-6) {
  vapply(seq_along(x), function(j) {
    xp <- x
    xm <- x
    xp[j] <- xp[j] + eps
    xm[j] <- xm[j] - eps
    (f(xp) - f(xm)) / (2 * eps)
  }, numeric(1))
}

test_that("PH cumhaz beta/gamma gradients agree with finite differences", {
  dat <- survival::veteran
  fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = dat, approach = "mle", init = 0)
  b <- fit$coefficients
  g <- fit$bp.param
  m <- length(g)
  t0 <- 50
  xrow <- c(karno = 60, "factor(celltype)adeno" = 0, "factor(celltype)smallcell" = 0, "factor(celltype)squamous" = 1)
  G <- sapply(seq_len(m), function(k) pbeta(t0 / fit$tau_b, k, m - k + 1))
  eta <- function(beta) as.numeric(exp(sum(xrow * beta)))
  H <- function(beta, gamma) as.numeric(sum(G * gamma) * eta(beta))

  num_b <- .fd_grad(function(beta) H(beta, g), b)
  num_g <- .fd_grad(function(gamma) H(b, gamma), g)

  exp_eta <- eta(b)
  ana_b <- as.numeric(sum(G * g) * exp_eta) * xrow
  ana_g <- G * exp_eta

  expect_equal(unname(num_b), unname(ana_b), tolerance = 1e-5)
  expect_equal(unname(num_g), unname(ana_g), tolerance = 1e-5)
})

test_that("PO cumhaz beta/gamma gradients agree with finite differences", {
  dat <- survival::veteran
  fit <- bppo(Surv(time, status) ~ karno + factor(celltype), data = dat, approach = "mle", init = 0)
  b <- fit$coefficients
  g <- fit$bp.param
  m <- length(g)
  t0 <- 50
  xrow <- c(karno = 60, "factor(celltype)adeno" = 0, "factor(celltype)smallcell" = 0, "factor(celltype)squamous" = 1)
  G <- sapply(seq_len(m), function(k) pbeta(t0 / fit$tau_b, k, m - k + 1))
  eta <- function(beta) as.numeric(exp(sum(xrow * beta)))
  H <- function(beta, gamma) {
    odds <- as.numeric(sum(G * gamma) * eta(beta))
    log(1 + odds)
  }

  num_b <- .fd_grad(function(beta) H(beta, g), b)
  num_g <- .fd_grad(function(gamma) H(b, gamma), g)

  exp_eta <- eta(b)
  odds <- as.numeric(sum(G * g) * exp_eta)
  ana_b <- as.numeric(odds / (1 + odds)) * xrow
  ana_g <- as.numeric(1 / (1 + odds)) * G * exp_eta

  expect_equal(unname(num_b), unname(ana_b), tolerance = 1e-5)
  expect_equal(unname(num_g), unname(ana_g), tolerance = 1e-5)
})

test_that("AFT cumhaz beta/gamma gradients agree with finite differences (endpoint-corrected tau)", {
  veteran2 <- survival::veteran[survival::veteran$prior == 0, ]
  fit <- bpaft(Surv(time, status) ~ karno, data = veteran2, approach = "mle", degree = 5, init = 0)
  b <- fit$coefficients
  g <- fit$bp.param
  m <- length(g)
  t0 <- 100
  xrow <- c(karno = 60)
  X_train <- model.matrix(fit)
  log_time_train <- log(fit$y[, 1])
  P <- fit$standata$P
  if (is.null(P)) {
    P <- pw.basis(degree = m)
  }

  cumhaz_at <- function(beta, gamma) {
    spsurv:::.spbp_aft_cumhaz_grad(
      tau_a = fit$tau_a,
      tau_b = fit$tau_b,
      t = t0,
      gamma = gamma,
      beta = beta,
      xrow = xrow,
      P = P,
      X_train = X_train,
      log_time_train = log_time_train
    )$cumhaz
  }

  num_b <- .fd_grad(function(beta) cumhaz_at(beta, g), b)
  num_g <- .fd_grad(function(gamma) cumhaz_at(b, gamma), g)

  ana <- spsurv:::.spbp_aft_cumhaz_grad(
    tau_a = fit$tau_a,
    tau_b = fit$tau_b,
    t = t0,
    gamma = g,
    beta = b,
    xrow = xrow,
    P = P,
    X_train = X_train,
    log_time_train = log_time_train
  )

  expect_equal(unname(num_b), unname(ana$grad_beta), tolerance = 1e-3)
  expect_equal(unname(num_g), unname(ana$grad_gamma), tolerance = 1e-3)
})

test_that("BPAFT production survfit cumhaz matches endpoint-corrected map and FD gradients", {
  veteran2 <- survival::veteran[survival::veteran$prior == 0, ]
  fit <- bpaft(Surv(time, status) ~ karno, data = veteran2, approach = "mle", degree = 5, init = 0)
  b <- fit$coefficients
  g <- fit$bp.param
  m <- length(g)
  t0 <- 75
  nd <- data.frame(karno = 60)
  X_train <- model.matrix(fit)
  log_time_train <- log(fit$y[, 1])
  P <- fit$standata$P
  if (is.null(P)) {
    P <- pw.basis(degree = m)
  }

  sf <- survfit(fit, newdata = nd, times = t0)
  xrow <- as.numeric(stats::model.matrix(fit$formula[-2], xlev = fit$xlevels, data = nd)[, -1, drop = TRUE])
  ana <- spsurv:::.spbp_aft_cumhaz_grad(
    tau_a = fit$tau_a,
    tau_b = fit$tau_b,
    t = t0,
    gamma = g,
    beta = b,
    xrow = xrow,
    P = P,
    X_train = X_train,
    log_time_train = log_time_train
  )

  expect_equal(unname(sf$cumhaz[1, 1]), unname(ana$cumhaz), tolerance = 1e-10)
  expect_equal(unname(sf$surv[1, 1]), exp(-ana$cumhaz), tolerance = 1e-10)

  cumhaz_at <- function(beta, gamma) {
    spsurv:::.spbp_aft_cumhaz_grad(
      tau_a = fit$tau_a,
      tau_b = fit$tau_b,
      t = t0,
      gamma = gamma,
      beta = beta,
      xrow = xrow,
      P = P,
      X_train = X_train,
      log_time_train = log_time_train
    )$cumhaz
  }
  num_b <- .fd_grad(function(beta) cumhaz_at(beta, g), b)
  num_g <- .fd_grad(function(gamma) cumhaz_at(b, gamma), g)
  expect_equal(unname(num_b), unname(ana$grad_beta), tolerance = 1e-3)
  expect_equal(unname(num_g), unname(ana$grad_gamma), tolerance = 1e-3)

  if (.spbp_gamma_information_diagnostics(fit)$stable) {
    V <- vcov(fit, bp.param = TRUE)
    grad <- c(ana$grad_beta, ana$grad_gamma)
    se_manual <- sqrt(as.numeric(t(grad) %*% V %*% grad))
    expect_equal(unname(sf$std.err[1, 1]), unname(se_manual), tolerance = 1e-8)
  }
})

test_that("PO survfit std.err follows delta method when gamma information is stable", {
  dat <- survival::veteran
  # Fixed low degree + init=0 keeps the gamma information finite and well-conditioned
  # across platforms (default degree ceiling(sqrt(n)) is often ill-conditioned for PO).
  fit <- bppo(
    Surv(time, status) ~ karno + factor(celltype),
    data = dat,
    approach = "mle",
    degree = 5L,
    init = 0
  )
  t0 <- 50
  gamma_diag <- .spbp_gamma_information_diagnostics(fit)
  expect_true(gamma_diag$available)
  expect_true(gamma_diag$stable)
  expect_true(is.finite(gamma_diag$kappa_gamma))

  sf <- survfit(fit, times = t0)
  q <- length(fit$coefficients)
  m <- length(fit$bp.param)
  X <- t(fit$means)
  G <- sapply(seq_len(m), function(k) pbeta(t0 / fit$tau_b, k, m - k + 1))
  exp_eta <- exp(as.numeric(X %*% fit$coefficients))
  odds <- sum(G * fit$bp.param) * exp_eta
  H <- log1p(odds)
  grad <- c(odds / (1 + odds) * as.numeric(X), G * exp_eta / (1 + odds))
  V <- vcov(fit, bp.param = TRUE)
  se_manual <- sqrt(as.numeric(t(grad) %*% V %*% grad))
  expect_equal(unname(sf$std.err[1, 1]), unname(se_manual), tolerance = 1e-8)
  expect_lt(se_manual / H, 1)
})

test_that("non-finite gamma information does not crash diagnostics or summary/survfit", {
  dat <- survival::veteran
  fit <- bppo(
    Surv(time, status) ~ karno + factor(celltype),
    data = dat,
    approach = "mle",
    degree = 5L,
    init = 0
  )
  q <- length(fit$coefficients)
  # Inject a non-finite entry into the gamma block of the stored Hessian.
  bad <- fit
  bad$hessian[q + 1L, q + 1L] <- Inf

  gamma_diag <- .spbp_gamma_information_diagnostics(bad)
  expect_false(gamma_diag$stable)
  expect_false(gamma_diag$available)
  expect_true(is.na(gamma_diag$rank_gamma))
  expect_true(is.na(gamma_diag$kappa_gamma))
  expect_match(gamma_diag$reason, "non-finite")

  expect_warning(s <- summary(bad), regexp = NULL)
  expect_true(is.matrix(s$coefficients))
  expect_true(all(is.finite(s$coefficients[, "se(coef)"])))

  expect_warning(sf <- survfit(bad, times = 50), "non-finite|Bernstein")
  expect_true(all(is.na(sf$std.err)))
  expect_true(all(is.finite(sf$surv)))
})

test_that("survfit point estimates satisfy survival identity", {
  dat <- survival::veteran
  fit <- bpph(Surv(time, status) ~ karno, data = dat, approach = "mle", init = 0)
  sf <- survfit(fit, times = seq(0, 120, by = 10))
  expect_equal(sf$surv, exp(-sf$cumhaz), tolerance = 1e-10)
})

test_that("AFT power-basis helper matches beta-mixture basis at fixed range", {
  dat <- survival::veteran
  fit <- bpaft(Surv(time, status) ~ karno, data = dat, approach = "mle", init = 0)
  p <- fit$standata$P
  u <- seq(0.05, 0.95, by = 0.05)
  y <- fit$tau_a + u * (fit$tau_b - fit$tau_a)
  tgrid <- exp(y)
  lp <- rep(0, length(tgrid))
  y <- log(tgrid) - lp
  m <- length(fit$bp.param)

  basis_pw <- spsurv:::.spbp_aft_basis(y, fit$tau_a, fit$tau_b, p)
  G_beta <- sapply(seq_len(m), function(k) pbeta(u, k, m - k + 1))
  g_beta <- sapply(seq_len(m), function(k) dbeta(u, k, m - k + 1)) / (fit$tau_b - fit$tau_a)

  expect_equal(unname(basis_pw$G), unname(G_beta), tolerance = 1e-10)
  expect_equal(unname(basis_pw$g), unname(g_beta), tolerance = 1e-10)
})
