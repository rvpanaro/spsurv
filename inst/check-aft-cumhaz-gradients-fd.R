#!/usr/bin/env Rscript
# BPAFT: analytic d(cumhaz)/d(beta, xi) vs FD through full tau(beta) recomputation.
# Run: Rscript inst/check-aft-cumhaz-gradients-fd.R

pkg_root <- Sys.getenv("SPSURV_PKG_ROOT", unset = normalizePath(getwd(), winslash = "/"))
if (file.exists(file.path(pkg_root, "DESCRIPTION"))) {
  suppressPackageStartupMessages(devtools::load_all(pkg_root, quiet = TRUE))
} else {
  library(spsurv)
}
library(survival)

.fd_grad <- function(f, x, eps = 1e-6) {
  vapply(seq_along(x), function(j) {
    xp <- x
    xm <- x
    xp[j] <- xp[j] + eps
    xm[j] <- xm[j] - eps
    (f(xp) - f(xm)) / (2 * eps)
  }, numeric(1))
}

.rel_err <- function(ana, num) {
  denom <- pmax(abs(num), abs(ana), 1e-12)
  abs(ana - num) / denom
}

.aft_cumhaz_full <- function(fit, t, xrow, beta, gamma, recompute_tau = TRUE) {
  m <- length(gamma)
  P <- fit$standata$P
  if (is.null(P)) {
    P <- pw.basis(degree = m)
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
    )$cumhaz
  } else {
    spsurv:::.spbp_aft_cumhaz_grad(fit$tau_a, fit$tau_b, t, gamma, beta, xrow, P)$cumhaz
  }
}

.aft_cumhaz_grad_ana <- function(fit, t, xrow, beta, gamma, recompute_tau = TRUE) {
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
    spsurv:::.spbp_aft_cumhaz_grad(fit$tau_a, fit$tau_b, t, gamma, beta, xrow, P)
  }
}

.check_case <- function(fit, t, xrow, x_label, eps = 1e-6, flag_tol = 1e-3) {
  beta <- fit$coefficients
  gamma <- fit$bp.param
  P <- fit$standata$P
  if (is.null(P)) {
    P <- pw.basis(degree = length(gamma))
  }

  ana <- .aft_cumhaz_grad_ana(fit, t, xrow, beta, gamma, recompute_tau = TRUE)

  num_beta_full <- .fd_grad(
    function(b) .aft_cumhaz_full(fit, t, xrow, b, gamma, recompute_tau = TRUE),
    beta,
    eps = eps
  )
  num_gamma_full <- .fd_grad(
    function(g) .aft_cumhaz_full(fit, t, xrow, beta, g, recompute_tau = TRUE),
    gamma,
    eps = eps
  )

  num_beta_frozen <- .fd_grad(
    function(b) .aft_cumhaz_full(fit, t, xrow, b, gamma, recompute_tau = FALSE),
    beta,
    eps = eps
  )

  beta_nms <- names(beta)
  gamma_nms <- names(gamma)
  if (is.null(beta_nms)) {
    beta_nms <- paste0("beta[", seq_along(beta), "]")
  }
  if (is.null(gamma_nms)) {
    gamma_nms <- paste0("xi[", seq_along(gamma), "]")
  }

  rows <- rbind(
    data.frame(
      parameter = beta_nms,
      block = "beta",
      analytic_code = unname(ana$grad_beta),
      numeric_full = unname(num_beta_full),
      numeric_frozen = unname(num_beta_frozen),
      rel_error_full = .rel_err(ana$grad_beta, num_beta_full),
      rel_error_frozen = .rel_err(ana$grad_beta, num_beta_frozen),
      stringsAsFactors = FALSE
    ),
    data.frame(
      parameter = gamma_nms,
      block = "xi",
      analytic_code = unname(ana$grad_gamma),
      numeric_full = unname(num_gamma_full),
      numeric_frozen = NA_real_,
      rel_error_full = .rel_err(ana$grad_gamma, num_gamma_full),
      rel_error_frozen = NA_real_,
      stringsAsFactors = FALSE
    )
  )
  rows$flagged_full <- rows$rel_error_full > flag_tol
  rows$flagged_frozen <- !is.na(rows$rel_error_frozen) & rows$rel_error_frozen > flag_tol

  attr(rows, "meta") <- list(
    t = t,
    x_label = x_label,
    cumhaz = ana$cumhaz,
    surv = exp(-ana$cumhaz),
    tau_a_fit = fit$tau_a,
    tau_b_fit = fit$tau_b
  )
  rows
}

print_case <- function(rows) {
  meta <- attr(rows, "meta")
  cat("\n--- t =", meta$t, "| x =", meta$x_label, "---\n")
  cat("cumhaz =", signif(meta$cumhaz, 6), " surv =", signif(meta$surv, 6), "\n")
  cat("frozen tau_a/tau_b:", signif(meta$tau_a_fit, 6), signif(meta$tau_b_fit, 6), "\n\n")
  print(rows, row.names = FALSE, digits = 8)
  cat("\nflagged (full FD):", sum(rows$flagged_full, na.rm = TRUE),
      "| flagged (frozen-tau FD):", sum(rows$flagged_frozen, na.rm = TRUE), "\n")
  invisible(rows)
}

veteran2 <- veteran[veteran$prior == 0, ]
veteran2$celltype <- factor(veteran2$celltype)

cat("=== veteran2 BPAFT fit (karno only, degree 5) ===\n")
fit <- bpaft(
  Surv(time, status) ~ karno,
  data = veteran2,
  approach = "mle",
  degree = 5,
  init = 0
)
gamma_diag <- spsurv:::.spbp_gamma_information_diagnostics(fit)
cat("stable:", gamma_diag$stable, "| kappa(C):", signif(gamma_diag$kappa_gamma, 4), "\n")
cat("tau_a:", fit$tau_a, " tau_b:", fit$tau_b, "\n\n")

times <- c(30, 75, 120)
x_profiles <- list(
  means = list(label = "training means", row = as.numeric(fit$means)),
  k60 = list(label = "karno=60", row = {
    nd <- data.frame(karno = 60)
    as.numeric(stats::model.matrix(fit$formula[-2], data = nd)[, -1, drop = TRUE])
  })
)

all_rows <- list()
for (t in times) {
  for (prof in x_profiles) {
    rows <- .check_case(fit, t, prof$row, prof$label)
    print_case(rows)
    all_rows[[length(all_rows) + 1L]] <- rows
  }
}

combined <- do.call(rbind, all_rows)
cat("\n========== SUMMARY ==========\n")
cat("Max rel_error (analytic vs FULL FD):", signif(max(combined$rel_error_full), 6), "\n")
cat("Max rel_error (analytic vs frozen-tau FD):", signif(max(combined$rel_error_frozen, na.rm = TRUE), 6), "\n")
cat("Flagged full:", sum(combined$flagged_full), "| frozen:", sum(combined$flagged_frozen, na.rm = TRUE), "\n")
if (any(combined$flagged_full)) {
  cat("\nFlagged (full):\n")
  print(subset(combined, flagged_full), row.names = FALSE)
}
