#---------------------------------------------
#' Covariance of the regression coefficients
#'
#' Uses block-wise inversion of the negative Hessian, with a clear split between
#' the regression coefficients (beta) and the Bernstein polynomial coefficients (gamma).
#'
#' @aliases vcov.spbp
#' @export
#' @param object an object of the class spbp
#' @param bp.param return Bernstein Polynomial variance.
#' @param ... arguments passed to parent method.
#' @return  the variance-covariance matrix associated with the regression coefficients.
#'

vcov.spbp <- function(object, bp.param = FALSE, ...) {
  if (object$call$approach != "mle") {
    warning("Not available, change approach to 'mle' instead.")
    return(NULL)
  }

  q <- length(object$coefficients)
  m <- length(object$bp.param)

  # Negative Hessian (observed information). Stan order: (alpha, beta, gamma)
  H_neg <- -object$hessian
  if (nrow(H_neg) == 1L + q + m) {
    # Drop alpha row/column so we work with (beta, gamma) only
    H_neg <- H_neg[-1L, -1L, drop = FALSE]
  }
  stopifnot(nrow(H_neg) == q + m)

  # Block-wise inversion: H = [A B; B' C], A = beta block, C = gamma block
  if (q == 0L) {
    # No regression coefficients (null model): only gamma block
    inv_C <- MASS::ginv(H_neg)
    V <- inv_C
  } else {
    A <- H_neg[1:q, 1:q, drop = FALSE]
    B <- H_neg[1:q, (q + 1L):(q + m), drop = FALSE]
    C <- H_neg[(q + 1L):(q + m), (q + 1L):(q + m), drop = FALSE]

    inv_C <- MASS::ginv(C)
    # Schur complement of C: S_A = A - B inv(C) B'
    Schur_A <- A - B %*% inv_C %*% t(B)
    V_bb <- MASS::ginv(Schur_A)
    V_bg <- -V_bb %*% B %*% inv_C
    V_gg <- inv_C + inv_C %*% t(B) %*% V_bb %*% B %*% inv_C

    V <- matrix(0, q + m, q + m)
    V[1:q, 1:q] <- V_bb
    V[1:q, (q + 1L):(q + m)] <- V_bg
    V[(q + 1L):(q + m), 1:q] <- t(V_bg)
    V[(q + 1L):(q + m), (q + 1L):(q + m)] <- V_gg
  }
  diag(V) <- abs(diag(V))

  # Parameters for Delta Method
  inv_exp_a <- 1 / exp(object$alpha)

  # Construct Jacobian J
  J <- matrix(0, q + m, q + m)

  # d(coef)/d(beta_std)
  diag(J)[1:q] <- 1 / object$sdv

  # d(gamma)/d(gamma_std)
  diag(J)[(q + 1):(q + m)] <- inv_exp_a

  # d(gamma)/d(beta_std): The cross-dependency
  # This accounts for how the intercept shift affects the gamma weights
  for (i in 1:m) {
    J[q + i, 1:q] <- -object$bp.param[i] * inv_exp_a * object$means
  }

  V_final <- J %*% V %*% t(J)

  # Labeling
  nms <- c(names(object$coefficients), names(object$bp.param))
  dimnames(V_final) <- list(nms, nms)

  if (bp.param) {
    return(V_final)
  } else {
    return(V_final[1:q, 1:q])
  }
}
