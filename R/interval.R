#---------------------------------------------
#' Generic S3 method confint
#' @aliases confint
#' @export
#' @param x a fitted model object
#' @param ... further arguments passed to or from other methods.
#' @return the estimated regression coefficients
#'
confint <- function(x, ...) UseMethod("confint")

#---------------------------------------------
#' Confidence intervals for the regression coefficients
#'
#' @aliases confint.spbp
#' @rdname confint-methods
#' @method confint spbp
#' @export
#' @export confint
#' @param x an object of the class x
#' @param level the confidence level required
#' @param ... further arguments passed to or from other methods
#' @return  100(1-alpha) confidence intervals for the regression coefficients
#'
#'
confint.spbp <- function(x, level = 0.95) {
  if (x$call$approach == "mle") {
    se <- sqrt(diag(as.matrix(vcov(x))))
    alpha <- 1 - level
    CI <- as.vector(coef(x)) + se %o% c(-qnorm(1 - alpha / 2), qnorm(1 - alpha / 2))

    labels <- round(100 * (c(alpha / 2, 1 - alpha / 2)), 1)
    colnames(CI) <- paste0(labels, "%")
    return(CI)
  } else {
    warning("not applicable, calling credint.spbp to get credible intervals instead")
    credint.spbp(x, prob = level)
  }
}
#---------------------------------------------
#' Generic S3 method credint
#' @aliases credint
#' @export
#' @param x a fitted model object
#' @param ... further arguments passed to or from other methods.
#' @return the estimated regression coefficients
#'
credint <- function(x, ...) UseMethod("credint")

#---------------------------------------------
#' Confidence intervals for the regression coefficients
#'
#' @aliases credint.spbp
#' @rdname credint-methods
#' @method credint spbp
#' @export
#' @export credint
#' @param x an object of the class x
#' @param prob the probability level required
#' @param ... further arguments passed to or from other methods
#' @return  100(1-alpha) confidence intervals for the regression coefficients
#'
#'
credint.spbp <- function(x, prob = 0.95, type = c("HPD", "Equal-Tailed")) {
  type <- match.arg(type)

  if (x$call$approach == "mle") {
    warning("not applicable, calling confint.spbp to get credible intervals instead")
    confint.spbp(x, level = prob)
  } else {
    if (type == "HPD") {
      coda::HPDinterval(coda::mcmc(x$posterior$beta))
    } else {
      apply(x$posterior$beta, 2, quantile, probs = c((1 - prob) / 2, 1 - (1 - prob) / 2)) %>% t()
    }
  }
}
