#' Bernstein Polynomial Based Regression Object Print
#'
#' @export
#' @param x an object of class spbp.
#' @param digits number of digits to display.
#' @param signif.stars see \code{\link{getOption}}.
#' @param bp.param print BP parameters.
#' @param ... further arguments passed to or from other methods.
#' @method print spbp
#' @return none
#' @importFrom stats pnorm

print.spbp <-
  function(x, bp.param = FALSE, digits = max(getOption("digits") - 4, 3),
           signif.stars = getOption("show.signif.stars"), ...) {
    savedig <- options(digits = digits)
    on.exit(options(savedig))

    if (!is.null(x$call)) {
      cat("Call:\n")
      dput(x$call)
      cat("\n")
    }

    coef <- coef(x)

    if (x$call$approach == "mle") {
      var <- vcov(x)

      ### Error handling ###
      # Null model
      coef2 <- coef[!(is.na(coef))] # non-missing coefs
      if (is.null(var)) stop("Input is not valid")

      if (!is.null(x$coefficients) & !bp.param) {
        se <- sqrt(diag(as.matrix(vcov(x))))

        Coefmat <- cbind(
          coef, exp(coef), se, coef / se,
          pchisq((coef / se)^2, 1, lower.tail = FALSE)
        )
        dimnames(Coefmat) <- list(names(coef), c(
          "coef", "exp(coef)",
          "se(coef)", "z", "Pr(>|z|)"
        ))
      } else {
        se <- sqrt(diag(as.matrix(vcov(x, TRUE))))

        Coefmat <- cbind(
          log(x$bp.param), x$bp.param, se, x$bp.param * (log(x$bp.param) - log(1e-19)) / se,
          pnorm(q = x$bp.param * (log(x$bp.param) - log(1e-19)) / se, lower.tail = FALSE)
        )

        dimnames(Coefmat) <- list(names(x$bp.param), c(
          "log(gamma)", "gamma",
          "se(log(gamma))", "z", "Pr(>z)"
        ))
      }

      cat("\n")
      printCoefmat(Coefmat,
        digits = digits,
        signif.stars = signif.stars, ...
      )

      if (!is.null(x$loglik)) {
        logtest <- -2 * (x$loglik[1] - x$loglik[2])

        cat(
          "\nLoglik(model)= ", x$loglik[2], "  Loglik(baseline only)= ", x$loglik[1], "\n"
        )
        if (!is.null(x$coefficients)) {
          cat(
            "        Chisq= ", logtest, " on ", ncol(x$features), " df, p= ",
            pchisq(logtest, ncol(x$features), lower.tail = FALSE), "\n"
          )
        }
      }

      if (!is.null(x$n)) {
        cat("n= ", paste0(x$n, ","), " number of events= ", x$nevent)
      }
    } else {
      coef2 <- coef[!(is.na(coef))] # non-missing coefs

      if (is.null(x$posterior)) stop("Input is not valid")


      if (!is.null(x$coefficients)) {
        Coefmat <- cbind(
          coef,
          apply(x$posterior$beta, 2, .mode),
          apply(x$posterior$beta, 2, median),
          colMeans(exp(x$posterior$beta)),
          apply(x$posterior$beta, 2, sd)
        )
        dimnames(Coefmat) <- list(names(coef), c(
          "mean(coef)",
          "mode(coef)",
          "median(coef)",
          "mean(exp(coef))", "sd(coef)"
        ))

        cat("\n")
        printCoefmat(Coefmat,
          digits = digits,
          signif.stars = signif.stars, ...
        )
      } else {
        Coefmat <- cbind(
          x$bp.param,
          apply(x$posterior$gamma, 2, .mode),
          apply(x$posterior$gamma, 2, median),
          colMeans(log(x$posterior$gamma)),
          apply(x$posterior$gamma, 2, sd)
        )
        dimnames(Coefmat) <- list(names(x$bp.param), c(
          "mean(bp)",
          "mode(bp)",
          "median(bp)",
          "mean(log(bp))", "sd(bp)"
        ))

        cat("\n")
        printCoefmat(Coefmat,
          digits = digits,
          signif.stars = signif.stars, ...
        )
      }
      cat("\n")

      if (!is.null(x$posterior$log_lik)) {
        cat(
          "\nDeviance criterion= ", DIC(x$posterior$log_lik)[1, 1], "  Watanabe\u2013Akaike criterion= ", WAIC(x$posterior$log_lik)[1, 1], "\nLog pseudo-marginal lik=", LPML(x$posterior$log_lik)[1], "\n"
        )
      }

      if (!is.null(x$n)) {
        cat("n= ", paste0(x$n, ","), " number of events= ", x$nevent)
      }
    }
  }
