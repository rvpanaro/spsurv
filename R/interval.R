#' Confidence intervals for the regression coefficients
#'
#' @aliases confint.spbp
#' @export
#' @param object a fitted model object.
#' @param parm a specification of which parameters are to be given confidence
#'   intervals: regression coefficient names and/or Bernstein baseline names
#'   (e.g. \code{"gamma[1]"}). If missing, all regression coefficients are used.
#' @param level the confidence level required.
#' @param ... further arguments passed to parent method
#' @return 100(1-alpha) confidence intervals for the requested parameters.
confint.spbp <- function(object, parm, level = 0.95, ...) {
  beta_nms <- names(coef(object))
  gamma_nms <- names(object$bp.param)
  has_covariates <- !is.null(object$features) && sum(object$features) != 0

  if (missing(parm)) {
    parm <- if (has_covariates) beta_nms else gamma_nms
  } else if (is.numeric(parm)) {
    all_nms <- c(beta_nms, gamma_nms)
    parm <- all_nms[parm]
  }

  beta_parm <- intersect(as.character(parm), beta_nms)
  gamma_parm <- intersect(as.character(parm), gamma_nms)

  if (!length(beta_parm) && !length(gamma_parm)) {
    stop("No matching parameters in 'parm'.", call. = FALSE)
  }

  alpha <- 1 - level
  crit <- c(-qnorm(1 - alpha / 2), qnorm(1 - alpha / 2))
  labels <- round(100 * (c(alpha / 2, 1 - alpha / 2)), 1)
  colnms <- paste0(labels, "%")

  pieces <- list()

  if (length(beta_parm)) {
    if (!has_covariates) {
      stop("No regression coefficients in this model.", call. = FALSE)
    }
    if (object$call$approach == "mle") {
      se <- sqrt(diag(as.matrix(vcov(object))))[beta_parm]
      est <- coef(object)[beta_parm]
      CI <- as.vector(est) + se %o% crit
      rownames(CI) <- beta_parm
      colnames(CI) <- colnms
      pieces$beta <- CI
    } else {
      warning("not applicable, calling credint.spbp to get credible intervals instead")
      return(credint.spbp(object, prob = level))
    }
  }

  if (length(gamma_parm)) {
    if (object$call$approach == "mle") {
      V <- vcov(object, bp.param = TRUE)
      se <- sqrt(diag(as.matrix(V)))[gamma_parm]
      est <- log(object$bp.param[gamma_parm])
      CI <- as.vector(est) + se %o% crit
      rownames(CI) <- gamma_parm
      colnames(CI) <- colnms
      pieces$gamma <- exp(CI)
    } else {
      warning("not applicable, calling credint.spbp to get credible intervals instead")
      return(credint.spbp(object, prob = level))
    }
  }

  if (length(pieces) == 1L) {
    return(pieces[[1L]])
  }

  CI <- do.call(rbind, pieces[intersect(c("beta", "gamma"), names(pieces))])
  CI[parm, , drop = FALSE]
}

#' Generic S3 method credint
#' @aliases credint
#' @export
#' @param x a fitted model object
#' @param ... further arguments passed to parent method
#' @return Generic function dispatching to \code{credint.spbp}.
#'
credint <- function(x, ...) UseMethod("credint")

#' Credible intervals for the regression coefficients
#'
#' @aliases credint.spbp
#' @export
#' @param x an object of the class x.
#' @param prob the probability level required.
#' @param type interval type.
#' @param ... further arguments passed to or from other methods.
#' @return 100(1-alpha) credible intervals for the regression coefficients
#' @importFrom stats quantile
credint.spbp <- function(x, prob = 0.95, type = c("HPD", "Equal-Tailed"), ...) {
  type <- match.arg(type)

  if (sum(x$features) != 0) {
    if (x$call$approach == "mle") {
      warning("not applicable, calling confint.spbp to get confidence intervals instead")
      confint.spbp(x, level = prob)
    } else {
      if (type == "HPD") {
        coda::HPDinterval(coda::mcmc(x$posterior$beta))
      } else {
        t(apply(x$posterior$beta, 2, quantile, probs = c((1 - prob) / 2, 1 - (1 - prob) / 2)))
      }
    }
  } else {
    if (x$call$approach == "mle") {
      warning("not applicable, calling confint.spbp to get confidence intervals instead")
      confint.spbp(x, level = prob)
    } else {
      if (type == "HPD") {
        coda::HPDinterval(coda::mcmc(x$posterior$gamma))
      } else {
        t(apply(x$posterior$gamma, 2, quantile, probs = c((1 - prob) / 2, 1 - (1 - prob) / 2)))
      }
    }
  }
}
