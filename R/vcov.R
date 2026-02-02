#---------------------------------------------
#' Generic S3 method vcov
#' @aliases vcov
#' @export
#' @param x a fitted model object
#' @param ... further arguments passed to or from other methods.
#' @return the variance-covariance matrix associated the regression coefficients.
#'
vcov <- function(x, ...) UseMethod("vcov")

#---------------------------------------------
#' Covariance of the regression coefficients
#'
#' @aliases vcov.spbp
#' @rdname vcov-methods
#' @method vcov spbp
#' @export
#' @export vcov
#' @param x an object of the class spbp
#' @param ... further arguments passed to or from other methods.
#' @return  the variance-covariance matrix associated with the regression coefficients.
#'

vcov.spbp <- function(x, bp.param.var = FALSE, ...) {
  if (x$call$approach != "mle") {
    warning("Not available, change approach to 'mle' instead.")
    return(NULL)
  }

  V <- MASS::ginv(-x$hessian)
  diag(V) <- abs(diag(V))

  q <- length(x$coefficients)
  m <- length(x$bp.param)

  # Parameters for Delta Method
  inv_exp_a <- 1 / exp(x$alpha)

  # Construct Jacobian J
  J <- matrix(0, q + m, q + m)

  # d(coef)/d(beta_std)
  diag(J)[1:q] <- 1 / x$sdv

  # d(gamma)/d(gamma_std)
  diag(J)[(q + 1):(q + m)] <- inv_exp_a

  # d(gamma)/d(beta_std): The cross-dependency
  # This accounts for how the intercept shift affects the gamma weights
  for (i in 1:m) {
    J[q + i, 1:q] <- -x$bp.param[i] * inv_exp_a * x$means
  }

  V_final <- J %*% V %*% t(J)

  # Labeling
  nms <- c(names(x$coefficients), names(x$bp.param))
  dimnames(V_final) <- list(nms, nms)

  if (bp.param.var) {
    return(V_final)
  } else {
    return(V_final[1:q, 1:q])
  }
}
