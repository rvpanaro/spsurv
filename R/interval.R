#' Confidence intervals for the regression coefficients
#'
#' @aliases confint.spbp
#' @export
#' @param object a fitted model object.
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the confidence level required.
#' @param ... further arguments passed to parent method
#' @return 100(1-alpha) confidence intervals for the regression coefficients
confint.spbp <- function(object, parm = names(coef(object)), level = 0.95, ...) {
  if (sum(object$features) != 0) {
    if (object$call$approach == "mle") {
      se <- sqrt(diag(as.matrix(vcov(object))))
      alpha <- 1 - level
      CI <- as.vector(coef(object)) + se %o% c(-qnorm(1 - alpha / 2), qnorm(1 - alpha / 2))

      labels <- round(100 * (c(alpha / 2, 1 - alpha / 2)), 1)
      colnames(CI) <- paste0(labels, "%")
      return(CI)
    } else {
      warning("not applicable, calling credint.spbp to get credible intervals instead")
      credint.spbp(object, prob = level)
    }
  } else {
    if (object$call$approach == "mle") {
      se <- sqrt(diag(as.matrix(vcov(object, bp.param = TRUE)))) / object$bp.param
      alpha <- 1 - level
      CI <- as.vector(log(object$bp.param)) + se %o% c(-qnorm(1 - alpha / 2), qnorm(1 - alpha / 2))

      labels <- round(100 * (c(alpha / 2, 1 - alpha / 2)), 1)
      colnames(CI) <- paste0(labels, "%")
      return(exp(CI))
    } else {
      warning("not applicable, calling credint.spbp to get credible intervals instead")
      credint.spbp(object, prob = level)
    }
  }
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
