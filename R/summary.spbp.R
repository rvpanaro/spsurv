#' Bernstein Polynomial Based Regression Object Summary
#'
#' @export
#' @param object an object of class spbp
#' @param interval interval coverage (confidence or credibility)
#' @param compact logical; if TRUE, print.summary methods show essential output only.
#' @param show_call logical; include the model call in printed summary.
#' @param show_intervals logical; include interval table in printed summary.
#' @param show_baseline logical; include Bernstein baseline (\code{gamma}) Wald
#'   table in printed summary (MLE only).
#' @param mle_test global test to print for MLE summaries: "lr" or "wald".
#' @param bayes_criterion global criterion to print for Bayesian summaries: "waic", "dic", or "lpml".
#' @param ... further arguments passed to or from other methods
#' @method summary spbp
#' @importFrom stats vcov
#' @return An object of class analogous to for e.g. 'summary.bppo.bayes'.

summary.spbp <- function(object, interval = 0.95,
                         compact = TRUE,
                         show_call = TRUE,
                         show_intervals = TRUE,
                         show_baseline = FALSE,
                         mle_test = c("lr", "wald"),
                         bayes_criterion = c("waic", "dic", "lpml"),
                         ...) {
  mle_test <- match.arg(mle_test)
  bayes_criterion <- match.arg(bayes_criterion)
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
    if (is.null(beta) || is.null(var)) stop("Input is not valid")

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
    output$coef_interval <- cbind(
      beta, beta - z * se, beta + z * se
    )
    dimnames(output$coef_interval) <- list(names(beta), c(
      "coef",
      paste("lower .", round(100 * interval, 2), sep = ""),
      paste("upper .", round(100 * interval, 2), sep = "")
    ))

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
    wald_stat <- coxph.wtest(vcov(object), as.numeric(unname(coef(object))))
    output$waldtest <- c(
      test = as.vector(round(wald_stat$test, 2)),
      df = df,
      pvalue = pchisq(as.vector(round(wald_stat$test, 2)), df,
        lower.tail = FALSE
      )
    )
    output$nparams <- .spbp_nparams(object)
    output$p <- ncol(vcov(object))
    output$controls <- list(
      compact = isTRUE(compact),
      show_call = isTRUE(show_call),
      show_intervals = isTRUE(show_intervals),
      show_baseline = isTRUE(show_baseline),
      mle_test = mle_test,
      bayes_criterion = bayes_criterion
    )

    if (isTRUE(show_baseline)) {
      output$baseline <- .spbp_baseline_mle_table(object)
    }

    class(output) <- switch(object$call$model,
      "po" = "summary.bppo.mle",
      "ph" = "summary.bpph.mle",
      "aft" = "summary.bpaft.mle"
    )
  } else {
    ### Error handling ###
    beta2 <- beta[!(is.na(beta))] # non-missing coefs
    if (is.null(beta) || is.null(object$posterior)) stop("Input is not valid")

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

    output$coef_interval <- cbind(
      beta,
      {
        hpd <- coda::HPDinterval(coda::mcmc(object$posterior$beta))
        cbind(hpd[, 1], hpd[, 2])
      }
    )
    dimnames(output$coef_interval) <- list(names(beta), c(
      "mean(coef)",
      paste("lower .", round(100 * interval, 2), "HPD", sep = ""),
      paste("upper .", round(100 * interval, 2), "HPD", sep = "")
    ))

    output$interval <- cbind(
      colMeans(exp(object$posterior$beta)),
      colMeans(exp(-object$posterior$beta)),
      {
        hpd <- coda::HPDinterval(coda::mcmc(exp(object$posterior$beta)))
        cbind(hpd[, 1], hpd[, 2])
      }
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
    output$controls <- list(
      compact = isTRUE(compact),
      show_call = isTRUE(show_call),
      show_intervals = isTRUE(show_intervals),
      show_baseline = isTRUE(show_baseline),
      mle_test = mle_test,
      bayes_criterion = bayes_criterion
    )

    class(output) <- switch(object$call$model,
      "po" = "summary.bppo.bayes",
      "ph" = "summary.bpph.bayes",
      "aft" = "summary.bpaft.bayes"
    )
  }
  output
}
