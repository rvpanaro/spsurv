# Finite-difference checks for BPPO R(t|x) gradients

library(spsurv)
library(survival)

context("BPPO R(t|x) gradient finite-difference check")

.fd_grad <- function(f, x, eps = 1e-6) {
  vapply(seq_along(x), function(j) {
    xp <- x
    xm <- x
    xp[j] <- xp[j] + eps
    xm[j] <- xm[j] - eps
    (f(xp) - f(xm)) / (2 * eps)
  }, numeric(1))
}

.po_R <- function(tau_b, t, gamma, beta, xrow) {
  .spbp_po_R_grad(tau_b, t, gamma, beta, xrow)$R
}

test_that("BPPO dR/dbeta and dR/dxi match central FD on veteran2", {
  veteran2 <- veteran[veteran$prior == 0, ]
  veteran2$celltype <- factor(veteran2$celltype)

  fit <- bppo(
    Surv(time, status) ~ karno + factor(celltype),
    data = veteran2,
    approach = "mle",
    degree = 5,
    init = 0
  )
  beta <- fit$coefficients
  gamma <- fit$bp.param

  times <- c(30, 75, 120)
  x_profiles <- list(
    as.numeric(fit$means),
    {
      nd <- data.frame(karno = 60, celltype = "squamous")
      as.numeric(stats::model.matrix(fit$formula[-2], xlev = fit$xlevels, data = nd)[, -1, drop = TRUE])
    }
  )

  max_rel <- 0
  for (t in times) {
    for (xrow in x_profiles) {
      ana <- .spbp_po_R_grad(fit$tau_b, t, gamma, beta, xrow)
      num_b <- .fd_grad(function(b) .po_R(fit$tau_b, t, gamma, b, xrow), beta)
      num_g <- .fd_grad(function(g) .po_R(fit$tau_b, t, g, beta, xrow), gamma)
      rel_b <- abs(ana$grad_beta - num_b) / pmax(abs(num_b), abs(ana$grad_beta), 1e-12)
      rel_g <- abs(ana$grad_gamma - num_g) / pmax(abs(num_g), abs(ana$grad_gamma), 1e-12)
      max_rel <- max(max_rel, rel_b, rel_g)
      expect_lt(max(rel_b, rel_g), 1e-3)
    }
  }
})
