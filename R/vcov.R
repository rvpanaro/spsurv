#' Covariance of the regression coefficients
#'
#' Uses block-wise inversion of the negative Hessian, with a clear split between
#' the regression coefficients (beta) and the Bernstein polynomial coefficients (gamma).
#'
#' @aliases vcov.spbp
#' @export
#' @param object an object of the class spbp
#' @param bp.param return Bernstein Polynomial variance.
#' @param mask_unstable_gamma when \code{TRUE} (default), set
#'   \eqn{\gamma}-related entries to \code{NA} if the Bernstein information block
#'   is ill-conditioned. Set \code{FALSE} internally when propagating survival
#'   uncertainty despite instability.
#' @param ... arguments passed to parent method.
#' @return  the variance-covariance matrix associated with the regression coefficients.
#'

vcov.spbp <- function(object, bp.param = FALSE, mask_unstable_gamma = TRUE, ...) {
  if (object$call$approach != "mle") {
    warning("Not available, change approach to 'mle' instead.")
    return(NULL)
  }

  q <- length(object$coefficients)
  m <- length(object$bp.param)
  gamma_diag <- .spbp_gamma_information_diagnostics(object)

  hx <- .spbp_hessian_neg_beta_gamma(object)
  H_neg <- hx$H_neg
  if (is.null(H_neg)) {
    warning("Hessian is not available; cannot compute vcov.")
    return(NULL)
  }

  # Block-wise inversion: H = [A B; B' C], A = beta block, C = gamma block
  gamma_unavailable <- isFALSE(gamma_diag$available)
  if (q == 0L) {
    # No regression coefficients (null model): only gamma block
    if (gamma_unavailable) {
      V <- matrix(NA_real_, m, m)
    } else {
      inv_C <- tryCatch(
        .spbp_sym_inv(H_neg),
        error = function(e) NULL
      )
      if (is.null(inv_C)) {
        warning(
          "Could not invert the gamma information block; returning NA variances.",
          call. = FALSE
        )
        V <- matrix(NA_real_, m, m)
      } else {
        V <- inv_C
      }
    }
  } else {
    A <- H_neg[1:q, 1:q, drop = FALSE]
    B <- H_neg[1:q, (q + 1L):(q + m), drop = FALSE]
    C <- H_neg[(q + 1L):(q + m), (q + 1L):(q + m), drop = FALSE]

    if (gamma_unavailable || any(!is.finite(C))) {
      # Beta block remains usable; gamma-related uncertainty is unavailable.
      V_bb <- tryCatch(
        .spbp_sym_inv(A),
        error = function(e) NULL
      )
      if (is.null(V_bb)) {
        warning(
          "Could not invert the beta information block; cannot compute vcov.",
          call. = FALSE
        )
        return(NULL)
      }
      V <- matrix(NA_real_, q + m, q + m)
      V[1:q, 1:q] <- V_bb
      warning(
        "Bernstein-polynomial (gamma) information matrix is unavailable or non-finite; ",
        "gamma-related variance entries are set to NA.",
        call. = FALSE
      )
    } else {
      inv_C <- tryCatch(
        .spbp_sym_inv(C),
        error = function(e) NULL
      )
      if (is.null(inv_C)) {
        warning(
          "Could not invert the gamma information block; gamma-related variance entries are set to NA.",
          call. = FALSE
        )
        V_bb <- tryCatch(
          .spbp_sym_inv(A),
          error = function(e) NULL
        )
        if (is.null(V_bb)) {
          warning(
            "Could not invert the beta information block; cannot compute vcov.",
            call. = FALSE
          )
          return(NULL)
        }
        V <- matrix(NA_real_, q + m, q + m)
        V[1:q, 1:q] <- V_bb
      } else {
        # Schur complement of C: S_A = A - B inv(C) B'
        Schur_A <- A - B %*% inv_C %*% t(B)
        V_bb <- .spbp_sym_inv(Schur_A)
        V_bg <- -V_bb %*% B %*% inv_C
        V_gg <- inv_C + inv_C %*% t(B) %*% V_bb %*% B %*% inv_C

        V <- matrix(0, q + m, q + m)
        V[1:q, 1:q] <- V_bb
        V[1:q, (q + 1L):(q + m)] <- V_bg
        V[(q + 1L):(q + m), 1:q] <- t(V_bg)
        V[(q + 1L):(q + m), (q + 1L):(q + m)] <- V_gg
      }
    }
  }
  if (all(is.finite(diag(V)))) {
    diag(V) <- abs(diag(V))
  } else {
    d <- diag(V)
    d[is.finite(d)] <- abs(d[is.finite(d)])
    diag(V) <- d
  }

  # Parameters for Delta Method
  inv_exp_a <- 1 / exp(object$alpha)

  # Construct Jacobian J
  J <- matrix(0, q + m, q + m)

  # d(coef)/d(beta_std)
  if (q > 0L) {
    diag(J)[1:q] <- 1 / object$sdv
  }

  # d(gamma)/d(gamma_std)
  if (m > 0L) {
    diag(J)[(q + 1):(q + m)] <- inv_exp_a
  }

  # d(gamma)/d(beta_std): covariate centering shifts the baseline (PH/PO only).
  # gamma_k = psi_k / exp(alpha), alpha = sum(beta_std_j * means_j / sdv_j)
  # => d(gamma_k)/d(beta_std_j) = -gamma_k * means_j / sdv_j
  # Skip the coupling when gamma information is unavailable so NA gamma blocks
  # do not contaminate finite beta variances under the delta-method map.
  if (!gamma_unavailable && object$call$model %in% c("ph", "po") && q > 0L) {
    for (i in seq_len(m)) {
      J[q + i, 1:q] <- -object$bp.param[i] * object$means / object$sdv
    }
  }

  if (gamma_unavailable && q > 0L) {
    # Map only the finite beta block; leave gamma-related entries as NA.
    J_bb <- J[1:q, 1:q, drop = FALSE]
    V_bb <- V[1:q, 1:q, drop = FALSE]
    V_final <- matrix(NA_real_, q + m, q + m)
    V_final[1:q, 1:q] <- J_bb %*% V_bb %*% t(J_bb)
  } else if (gamma_unavailable && q == 0L) {
    V_final <- matrix(NA_real_, m, m)
  } else {
    V_final <- J %*% V %*% t(J)
  }

  # Labeling
  nms <- c(names(object$coefficients), names(object$bp.param))
  dimnames(V_final) <- list(nms, nms)
  attr(V_final, "gamma_information_stable") <- gamma_diag$stable
  attr(V_final, "gamma_information_kappa") <- gamma_diag$kappa_gamma
  attr(V_final, "gamma_information_available") <- isTRUE(gamma_diag$available)

  if ((!gamma_diag$stable || gamma_unavailable) && isTRUE(mask_unstable_gamma)) {
    if (bp.param) {
      if (q > 0L) {
        V_final[(q + 1L):(q + m), ] <- NA_real_
        V_final[, (q + 1L):(q + m)] <- NA_real_
      } else {
        V_final[] <- NA_real_
      }
    }
  }

  if (bp.param) {
    V_final
  } else {
    V_final[1:q, 1:q, drop = FALSE]
  }
}
