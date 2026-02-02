#' Bernstein Polynomial Based Regression Object Summary Bayes
#'
#' @export
#' @param x a summary.spbp.bayes object
#' @param digits number of digits to display.
#' @param ... further arguments passed to or from other methods
#' @method print summary.spbp.bayes
#' @return none

print.summary.spbp.bayes <- ## summary printings
  function(x, digits = max(getOption("digits") - 4, 3), ...) {
    if (!is.null(x$call)) {
      cat("Call:\n")
      dput(x$call)
      cat("\n")
    }

    savedig <- options(digits = digits)
    on.exit(options(savedig))

    cat("  n=", x$n)
    if (!is.null(x$nevent)) {
      cat(", number of events=", x$nevent, "\n")
    } else {
      cat("\n")
    }

    if (nrow(x$coef) == 0) { # Null model
      cat("Null model\n")
      return()
    }

    if (!is.null(x$coefficients)) {
      cat("\n")
      printCoefmat(x$coefficients,
        digits = digits,
        signif.stars = signif.stars, ...
      )
    }
    if (!is.null(x$interval)) {
      cat("---\n")
      print(x$interval)
    }
    cat("\n")

    pdig <- max(1, getOption("digits") - 4) # default it too high IMO

    cat(
      "DIC= ", x$dic, "  WAIC= ", x$waic, " LPML=", x$lpml, "\n"
    )

    invisible()
  }
