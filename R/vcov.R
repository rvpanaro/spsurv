#---------------------------------------------
#' Covariance of the regression coefficients
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

  V <- MASS::ginv(-object$hessian)
  diag(V) <- abs(diag(V))

  q <- length(object$coefficients)
  m <- length(object$bp.param)

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
