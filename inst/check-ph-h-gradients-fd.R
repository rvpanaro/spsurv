#!/usr/bin/env Rscript
# Check BPPH analytic dH/d(beta), dH/d(gamma) vs central finite differences.
# Run from package root: Rscript inst/check-ph-h-gradients-fd.R

pkg_root <- Sys.getenv("SPSURV_PKG_ROOT", unset = normalizePath(getwd(), winslash = "/"))
if (file.exists(file.path(pkg_root, "DESCRIPTION"))) {
  suppressPackageStartupMessages(devtools::load_all(pkg_root, quiet = TRUE))
} else {
  library(spsurv)
}

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

.ph_H <- function(tau_b, t, gamma, beta, xrow) {
  spsurv:::.spbp_ph_cumhaz_grad(tau_b, t, gamma, beta, xrow)$H
}

.simulate_bpph <- function(n = 250, degree = 5, tau_b = 150, seed = 20260701L) {
  set.seed(seed)
  m <- degree
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

  list(
    data = data.frame(time = time, status = status, x1 = x1, x2 = x2),
    beta_true = beta_true,
    gamma_true = gamma_true,
    tau_b = tau_b,
    degree = degree
  )
}

.check_case <- function(fit, t, xrow, x_label, eps = 1e-6, flag_tol = 1e-3) {
  beta <- fit$coefficients
  gamma <- fit$bp.param
  tau_b <- fit$tau_b

  ana <- spsurv:::.spbp_ph_cumhaz_grad(tau_b, t, gamma, beta, xrow)
  H <- ana$H

  num_beta <- if (length(beta)) {
    .fd_grad(function(b) .ph_H(tau_b, t, gamma, b, xrow), beta, eps = eps)
  } else {
    numeric(0)
  }
  num_gamma <- .fd_grad(function(g) .ph_H(tau_b, t, g, beta, xrow), gamma, eps = eps)

  beta_nms <- if (length(beta)) names(beta) else character(0)
  gamma_nms <- names(gamma)
  if (is.null(beta_nms)) {
    beta_nms <- paste0("beta[", seq_along(beta), "]")
  }
  if (is.null(gamma_nms)) {
    gamma_nms <- paste0("gamma[", seq_along(gamma), "]")
  }

  rows <- rbind(
    if (length(beta)) {
      data.frame(
        parameter = beta_nms,
        block = "beta",
        analytic = unname(ana$grad_beta),
        numeric = unname(num_beta),
        rel_error = .rel_err(ana$grad_beta, num_beta),
        flagged = .rel_err(ana$grad_beta, num_beta) > flag_tol,
        stringsAsFactors = FALSE
      )
    },
    data.frame(
      parameter = gamma_nms,
      block = "gamma",
      analytic = unname(ana$grad_gamma),
      numeric = unname(num_gamma),
      rel_error = .rel_err(ana$grad_gamma, num_gamma),
      flagged = .rel_err(ana$grad_gamma, num_gamma) > flag_tol,
      stringsAsFactors = FALSE
    )
  )

  attr(rows, "meta") <- list(
    t = t,
    x_label = x_label,
    H = H,
    xrow = xrow
  )
  rows
}

print_case <- function(rows) {
  meta <- attr(rows, "meta")
  cat("\n--- t =", meta$t, "| x =", meta$x_label, "| H(t|x) =", signif(meta$H, 6), "---\n")
  cat("xrow:", paste(format(meta$xrow, digits = 4), collapse = ", "), "\n\n")
  print(rows, row.names = FALSE, digits = 8)
  n_flag <- sum(rows$flagged)
  cat(if (n_flag) {
    paste0("\nFLAGGED: ", n_flag, " parameter(s) with rel_error > 1e-3\n")
  } else {
    "\nNo parameters flagged (all rel_error <= 1e-3).\n"
  })
  invisible(rows)
}

sim <- .simulate_bpph()
cat("Simulated BPPH data: n =", nrow(sim$data),
    "| degree =", sim$degree,
    "| tau_b =", sim$tau_b, "\n")
cat("True beta:", paste(names(sim$beta_true), signif(sim$beta_true, 4), sep = "=", collapse = ", "), "\n")
cat("True gamma:", paste(signif(sim$gamma_true, 4), collapse = ", "), "\n\n")

fit <- bpph(
  Surv(time, status) ~ x1 + x2,
  data = sim$data,
  approach = "mle",
  degree = sim$degree,
  init = 0
)

cat("Fitted beta:", paste(names(fit$coefficients), signif(fit$coefficients, 4), sep = "=", collapse = ", "), "\n")
cat("Fitted gamma:", paste(signif(fit$bp.param, 4), collapse = ", "), "\n")
cat("Fitted tau_b:", signif(fit$tau_b, 6), "\n")
gamma_diag <- spsurv:::.spbp_gamma_information_diagnostics(fit)
cat("Gamma information stable:", gamma_diag$stable,
    "| kappa =", signif(gamma_diag$kappa_gamma, 4), "\n")

times <- c(30, 75, 120)
p <- length(fit$coefficients)
x_profiles <- list(
  means = list(
    label = "training means (survfit default)",
    row = as.numeric(fit$means)
  ),
  low_x1 = list(
    label = "x1=45, x2=0",
    row = {
      nd <- data.frame(x1 = 45, x2 = 0)
      as.numeric(stats::model.matrix(fit$formula[-2], xlev = fit$xlevels, data = nd)[, -1, drop = TRUE])
    }
  )
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
cat("Cases:", length(all_rows), "(", length(times), "times x", length(x_profiles), "covariate vectors )\n")
cat("Total parameters checked:", nrow(combined), "\n")
cat("Max rel_error:", signif(max(combined$rel_error), 6), "\n")
cat("Flagged count:", sum(combined$flagged), "\n")
if (any(combined$flagged)) {
  cat("\nFlagged entries:\n")
  print(subset(combined, flagged), row.names = FALSE, digits = 8)
}
