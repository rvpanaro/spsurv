# BPAFT cumhaz gradients: analytic vs FD

library(spsurv)
library(survival)

context("BPAFT cumhaz gradient finite-difference check")

.fd_grad <- function(f, x, eps = 1e-6) {
  vapply(seq_along(x), function(j) {
    xp <- x
    xm <- x
    xp[j] <- xp[j] + eps
    xm[j] <- xm[j] - eps
    (f(xp) - f(xm)) / (2 * eps)
  }, numeric(1))
}

.aft_cumhaz_at <- function(fit, t, xrow, beta, gamma, recompute_tau = FALSE) {
  m <- length(gamma)
  P <- fit$standata$P
  if (is.null(P)) {
    P <- pw.basis(degree = m)
  }
  if (recompute_tau) {
    X_train <- model.matrix(fit)
    log_time_train <- log(fit$y[, 1])
    spsurv:::.spbp_aft_cumhaz_grad(
      tau_a = fit$tau_a,
      tau_b = fit$tau_b,
      t = t,
      gamma = gamma,
      beta = beta,
      xrow = xrow,
      P = P,
      X_train = X_train,
      log_time_train = log_time_train
    )$cumhaz
  } else {
    spsurv:::.spbp_aft_cumhaz_grad(
      fit$tau_a, fit$tau_b, t, gamma, beta, xrow, P
    )$cumhaz
  }
}

.aft_cumhaz_grad <- function(fit, t, xrow, beta, gamma, recompute_tau = FALSE) {
  P <- fit$standata$P
  if (is.null(P)) {
    P <- pw.basis(degree = length(gamma))
  }
  if (recompute_tau) {
    spsurv:::.spbp_aft_cumhaz_grad(
      tau_a = fit$tau_a,
      tau_b = fit$tau_b,
      t = t,
      gamma = gamma,
      beta = beta,
      xrow = xrow,
      P = P,
      X_train = model.matrix(fit),
      log_time_train = log(fit$y[, 1])
    )
  } else {
    spsurv:::.spbp_aft_cumhaz_grad(
      fit$tau_a, fit$tau_b, t, gamma, beta, xrow, P
    )
  }
}

test_that("BPAFT analytic grad matches FD when tau_a/tau_b are frozen at fit values", {
  veteran2 <- veteran[veteran$prior == 0, ]
  fit <- bpaft(
    Surv(time, status) ~ karno,
    data = veteran2,
    approach = "mle",
    degree = 5,
    init = 0
  )
  beta <- fit$coefficients
  gamma <- fit$bp.param

  for (t in c(30, 75, 120)) {
    for (xrow in list(as.numeric(fit$means), c(karno = 60))) {
      ana <- .aft_cumhaz_grad(fit, t, xrow, beta, gamma, recompute_tau = FALSE)
      num_b <- .fd_grad(function(b) .aft_cumhaz_at(fit, t, xrow, b, gamma), beta)
      num_g <- .fd_grad(function(g) .aft_cumhaz_at(fit, t, xrow, beta, g), gamma)
      expect_equal(unname(ana$grad_beta), unname(num_b), tolerance = 1e-5)
      expect_equal(unname(ana$grad_gamma), unname(num_g), tolerance = 1e-5)
    }
  }
})

test_that("BPAFT analytic grad matches full tau(beta) recomputation FD for beta", {
  veteran2 <- veteran[veteran$prior == 0, ]
  fit <- bpaft(Surv(time, status) ~ karno, data = veteran2, approach = "mle", degree = 5, init = 0)
  beta <- fit$coefficients
  gamma <- fit$bp.param
  xrow <- c(karno = 60)

  for (t in c(30, 75, 120)) {
    ana <- .aft_cumhaz_grad(fit, t, xrow, beta, gamma, recompute_tau = TRUE)
    num_b <- .fd_grad(
      function(b) .aft_cumhaz_at(fit, t, xrow, b, gamma, recompute_tau = TRUE),
      beta
    )
    num_g <- .fd_grad(
      function(g) .aft_cumhaz_at(fit, t, xrow, beta, g, recompute_tau = TRUE),
      gamma
    )
    expect_equal(unname(ana$grad_beta), unname(num_b), tolerance = 1e-3)
    expect_equal(unname(ana$grad_gamma), unname(num_g), tolerance = 1e-3)
  }
})
