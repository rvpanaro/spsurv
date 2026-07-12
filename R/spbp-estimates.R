#' Parameter estimates for fitted spbp models
#'
#' @param object A fitted \code{"spbp"} object.
#' @param ... Not used.
#' @return Named numeric vector of regression and Bernstein-polynomial parameters.
#' @export
estimates <- function(object, ...) {
  UseMethod("estimates")
}

#' @export
#' @method estimates spbp
estimates.spbp <- function(object, ...) {
  if (!inherits(object, "spbp")) {
    stop("Object must inherit from class 'spbp'.", call. = FALSE)
  }

  if (identical(object$call$approach, "mle")) {
    beta <- coef(object)
    gamma <- object$bp.param
    if (is.null(beta)) {
      beta <- numeric(0)
    }
    if (is.null(gamma)) {
      gamma <- numeric(0)
    }
    if (is.null(names(gamma)) && length(gamma) > 0L) {
      names(gamma) <- paste0("gamma", seq_along(gamma))
    }
    return(c(beta, gamma))
  }

  if (is.null(object$posterior)) {
    return(NULL)
  }

  out <- numeric(0)
  if (!is.null(object$posterior$beta)) {
    beta <- apply(object$posterior$beta, 2, mean)
    out <- c(out, beta)
  }
  if (!is.null(object$posterior$gamma)) {
    gamma <- apply(object$posterior$gamma, 2, mean)
    if (is.null(names(gamma))) {
      names(gamma) <- paste0("gamma", seq_along(gamma))
    }
    out <- c(out, gamma)
  }
  out
}

#' Standard errors for fitted spbp model parameters
#'
#' @param object A fitted \code{"spbp"} object.
#' @param ... Not used.
#' @return Named numeric vector of standard errors aligned with \code{\link{estimates}}.
#' @export
se <- function(object, ...) {
  UseMethod("se")
}

#' @export
#' @method se spbp
se.spbp <- function(object, ...) {
  if (!inherits(object, "spbp")) {
    stop("Object must inherit from class 'spbp'.", call. = FALSE)
  }

  if (identical(object$call$approach, "mle")) {
    est <- estimates(object)
    v <- vcov(object, bp.param = TRUE)
    if (is.null(v)) {
      return(rep(NA_real_, length(est)))
    }
    se_vals <- sqrt(diag(as.matrix(v)))
    names(se_vals) <- rownames(v)
    se_vals <- se_vals[names(est)]
    return(se_vals)
  }

  if (is.null(object$posterior)) {
    return(NULL)
  }

  out <- numeric(0)
  if (!is.null(object$posterior$beta)) {
    beta_se <- apply(object$posterior$beta, 2, stats::sd)
    out <- c(out, beta_se)
  }
  if (!is.null(object$posterior$gamma)) {
    gamma_se <- apply(object$posterior$gamma, 2, stats::sd)
    if (is.null(names(gamma_se))) {
      names(gamma_se) <- paste0("gamma", seq_along(gamma_se))
    }
    out <- c(out, gamma_se)
  }
  out
}
