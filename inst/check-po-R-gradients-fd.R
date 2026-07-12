#!/usr/bin/env Rscript
# Check BPPO analytic dR/d(beta), dR/d(xi) vs central finite differences.
# R(t|x) = odds in R/survfit.R (PO branch). Run: Rscript inst/check-po-R-gradients-fd.R

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

.po_R <- function(tau_b, t, gamma, beta, xrow) {
  spsurv:::.spbp_po_R_grad(tau_b, t, gamma, beta, xrow)$R
}

.check_case <- function(fit, t, xrow, x_label, eps = 1e-6, flag_tol = 1e-3) {
  beta <- fit$coefficients
  gamma <- fit$bp.param
  tau_b <- fit$tau_b

  ana <- spsurv:::.spbp_po_R_grad(tau_b, t, gamma, beta, xrow)
  R <- ana$R

  num_beta <- if (length(beta)) {
    .fd_grad(function(b) .po_R(tau_b, t, gamma, b, xrow), beta, eps = eps)
  } else {
    numeric(0)
  }
  num_gamma <- .fd_grad(function(g) .po_R(tau_b, t, g, beta, xrow), gamma, eps = eps)

  beta_nms <- if (length(beta)) names(beta) else character(0)
  gamma_nms <- names(gamma)
  if (is.null(beta_nms)) {
    beta_nms <- paste0("beta[", seq_along(beta), "]")
  }
  if (is.null(gamma_nms)) {
    gamma_nms <- paste0("xi[", seq_along(gamma), "]")
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
      block = "xi",
      analytic = unname(ana$grad_gamma),
      numeric = unname(num_gamma),
      rel_error = .rel_err(ana$grad_gamma, num_gamma),
      flagged = .rel_err(ana$grad_gamma, num_gamma) > flag_tol,
      stringsAsFactors = FALSE
    )
  )

  attr(rows, "meta") <- list(t = t, x_label = x_label, R = R, xrow = xrow)
  rows
}

print_case <- function(rows) {
  meta <- attr(rows, "meta")
  cat("\n--- t =", meta$t, "| x =", meta$x_label, "| R(t|x) =", signif(meta$R, 6), "---\n")
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

veteran2 <- veteran[veteran$prior == 0, ]
veteran2$celltype <- factor(veteran2$celltype)

cat("=== veteran2 BPPO fit ===\n")
fit <- bppo(
  Surv(time, status) ~ karno + factor(celltype),
  data = veteran2,
  approach = "mle",
  degree = 5,
  init = 0
)
gamma_diag <- spsurv:::.spbp_gamma_information_diagnostics(fit)
cat("stable:", gamma_diag$stable, "| kappa(C):", signif(gamma_diag$kappa_gamma, 4), "\n")
cat("Fitted beta:", paste(names(fit$coefficients), signif(fit$coefficients, 4), sep = "=", collapse = ", "), "\n")
cat("Fitted xi (bp.param):", paste(signif(fit$bp.param, 4), collapse = ", "), "\n\n")

times <- c(30, 75, 120)
x_profiles <- list(
  means = list(
    label = "training means (survfit default)",
    row = as.numeric(fit$means)
  ),
  karno60_squamous = list(
    label = "karno=60, celltype=squamous",
    row = {
      nd <- data.frame(karno = 60, celltype = "squamous")
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
cat("Cases:", length(all_rows), "\n")
cat("Total parameters checked:", nrow(combined), "\n")
cat("Max rel_error:", signif(max(combined$rel_error), 6), "\n")
cat("Flagged count:", sum(combined$flagged), "\n")
if (any(combined$flagged)) {
  cat("\nFlagged entries:\n")
  print(subset(combined, flagged), row.names = FALSE, digits = 8)
}
