#' Bernstein Polynomial Based Regression Object Summary
#'
#' @export
#' @param object an object of class spbp
#' @param interval interval coverage (confidence or credibility)
#' @param ... further arguments passed to or from other methods
#' @method summary spbp
#' @importFrom stats vcov
#' @return An object of class analogous to for e.g. 'summary.bppo.bayes'.

summary.spbp <- function(object, interval = 0.95, ...) {
  # Null model
  if (is.null(object$coefficients)) {
    return(object)
  }
  beta <- object$coefficients

  ## mle approach
  if (object$call$approach == "mle") {
    var <- as.matrix(vcov(object))

    ### Error handling ###
    beta2 <- beta[!(is.na(beta))] # non-missing coefs
    if (is.null(beta) | is.null(var)) stop("Input is not valid")

    se <- sqrt(diag(var))

    output <- list(
      call = object$call,
      return_code = object$return_code,
      n = object$n,
      loglik = object$loglik
    )

    if (!is.null(object$nevent)) {
      output$nevent <- object$nevent
    }

    output$coefficients <- cbind(
      beta, exp(beta), se, beta / se,
      pchisq((beta / se)^2, 1, lower.tail = FALSE)
    )
    dimnames(output$coefficients) <- list(names(beta), c(
      "coef", "exp(coef)",
      "se(coef)", "z", "Pr(>|z|)"
    ))

    z <- qnorm((1 + interval) / 2, 0, 1)
    output$interval <- cbind(
      exp(beta), exp(-beta), exp(beta - z * se),
      exp(beta + z * se)
    )
    dimnames(output$interval) <- list(names(beta), c(
      "exp(coef)", "exp(-coef)",
      paste("lower .", round(100 * interval, 2), sep = ""),
      paste("upper .", round(100 * interval, 2), sep = "")
    ))

    df <- length(beta2)
    logtest <- -2 * (object$loglik[1] - object$loglik[2])
    output$logtest <- c(
      test = logtest,
      df = df,
      pvalue = pchisq(logtest, df, lower.tail = FALSE)
    )
    output$rsq <- c(
      rsq = 1 - exp(-logtest / object$n),
      maxrsq = 1 - exp(2 * object$loglik[1] / object$n)
    )
    output$waldtest <- c(
      test = as.vector(round(coxph.wtest(vcov(object), coef(object))$test, 2)),
      df = df,
      pvalue = pchisq(as.vector(round(coxph.wtest(vcov(object), coef(object))$test, 2)), df,
        lower.tail = FALSE
      )
    )
    output$p <- ncol(vcov(object))

    class(output) <- switch(object$call$model,
      "po" = "summary.bppo.mle",
      "ph" = "summary.bpph.mle",
      "aft" = "summary.bpaft.mle"
    )
  } else {
    ### Error handling ###
    beta2 <- beta[!(is.na(beta))] # non-missing coefs
    if (is.null(beta) | is.null(object$posterior)) stop("Input is not valid")

    output <- list(
      call = object$call,
      n = object$n,
      loglik = object$loglik
    )

    if (!is.null(object$nevent)) {
      output$nevent <- object$nevent
    }

    output$coefficients <- cbind(
      beta,
      apply(object$posterior$beta, 2, median),
      colMeans(exp(object$posterior$beta)),
      apply(object$posterior$beta, 2, sd)
    )
    dimnames(output$coefficients) <- list(names(beta), c(
      "mean(coef)",
      "median(coef)", "mean(exp(coef))", "sd(coef)"
    ))

    output$interval <- cbind(
      colMeans(exp(object$posterior$beta)),
      colMeans(exp(-object$posterior$beta)),
      HPDinterval(coda::mcmc(exp(object$posterior$beta)))[, 1],
      HPDinterval(coda::mcmc(exp(object$posterior$beta)))[, 2]
    )
    dimnames(output$interval) <- list(names(beta), c(
      "mean(exp(coef))", "mean(exp(-coef))",
      paste("lower .", round(100 * interval, 2), "HPD", sep = ""),
      paste("upper .", round(100 * interval, 2), "HPD", sep = "")
    ))

    if (!is.null(object$posterior$log_lik)) {
      output$dic <- DIC(object$posterior$log_lik)[1, 1]
      output$waic <- WAIC(object$posterior$log_lik)[1, 1]
      output$lpml <- LPML(object$posterior$log_lik)[1]
    }

    class(output) <- switch(object$call$model,
      "po" = "summary.bppo.bayes",
      "ph" = "summary.bpph.bayes",
      "aft" = "summary.bpaft.bayes"
    )
  }
  return(output)
}
