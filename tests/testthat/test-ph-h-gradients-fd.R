# Analytic vs finite-difference gradients of H(t|x) for BPPH (PH model).
# Detailed numeric output: Rscript inst/check-ph-h-gradients-fd.R

library(spsurv)

context("BPPH H(t|x) gradient finite-difference check")

.fd_grad <- function(f, x, eps = 1e-6) {
  vapply(seq_along(x), function(j) {
    xp <- x
    xm <- x
    xp[j] <- xp[j] + eps
    xm[j] <- xm[j] - eps
    (f(xp) - f(xm)) / (2 * eps)
  }, numeric(1))
}

.ph_H <- function(tau_b, t, gamma, beta, xrow) {
  .spbp_ph_cumhaz_grad(tau_b, t, gamma, beta, xrow)$H
}

test_that("BPPH dH/dbeta and dH/dgamma match central FD at multiple (t, x)", {
  set.seed(20260701L)
  n <- 250
  m <- 5
  tau_b <- 150
  beta_true <- c(x1 = 0.025, x2 = -0.35)
  gamma_true <- seq(0.15, 0.85, length.out = m)

  H0 <- function(t, gamma) {
    G <- vapply(seq_len(m), function(k) stats::pbeta(t / tau_b, k, m - k + 1), numeric(1))
    sum(G * gamma)
  }
  S <- function(t, x, beta, gamma) {
    exp(-H0(t, gamma) * exp(sum(x * beta)))
  }
  gen_time <- function(x, beta, gamma) {
    u <- stats::runif(1)
    stats::uniroot(function(t) S(t, x, beta, gamma) - u, interval = c(1e-8, tau_b * 1.5))$root
  }

  x1 <- stats::rnorm(n, 60, 8)
  x2 <- stats::rbinom(n, 1, 0.45)
  cens <- stats::rexp(n, rate = 1 / (tau_b * 0.85))
  time <- numeric(n)
  status <- integer(n)
  for (i in seq_len(n)) {
    xrow <- c(x1[i], x2[i])
    ev <- gen_time(xrow, beta_true, gamma_true)
    time[i] <- min(ev, cens[i])
    status[i] <- as.integer(ev <= cens[i])
  }
  dat <- data.frame(time = time, status = status, x1 = x1, x2 = x2)

  fit <- bpph(Surv(time, status) ~ x1 + x2, data = dat, approach = "mle", degree = m, init = 0)
  beta <- fit$coefficients
  gamma <- fit$bp.param

  times <- c(30, 75, 120)
  x_profiles <- list(
    as.numeric(fit$means),
    {
      nd <- data.frame(x1 = 45, x2 = 0)
      as.numeric(stats::model.matrix(fit$formula[-2], xlev = fit$xlevels, data = nd)[, -1, drop = TRUE])
    }
  )

  max_rel <- 0
  for (t in times) {
    for (xrow in x_profiles) {
      ana <- .spbp_ph_cumhaz_grad(fit$tau_b, t, gamma, beta, xrow)
      num_b <- .fd_grad(function(b) .ph_H(fit$tau_b, t, gamma, b, xrow), beta)
      num_g <- .fd_grad(function(g) .ph_H(fit$tau_b, t, g, beta, xrow), gamma)
      rel_b <- abs(ana$grad_beta - num_b) / pmax(abs(num_b), abs(ana$grad_beta), 1e-12)
      rel_g <- abs(ana$grad_gamma - num_g) / pmax(abs(num_g), abs(ana$grad_gamma), 1e-12)
      max_rel <- max(max_rel, rel_b, rel_g)
      expect_lt(max(rel_b, rel_g), 1e-3)
    }
  }
})
