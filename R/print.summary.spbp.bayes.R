#' Bernstein Polynomial based Regression Object
#'
#' @export
#' @method print summary.spbp.bayes
#'

print.summary.spbp.bayes <- ## summary printings
  function(x, digits = max(getOption('digits')-4, 3)){
    if (!is.null(x$call)) {
      cat("Call:\n")
      dput(x$call)
      cat("\n")
    }
    savedig <- options(digits = digits)
    on.exit(options(savedig))

    cat("  n=", x$n)
    if (!is.null(x$nevent)) cat(", number of events=", x$nevent, "\n")
    else cat("\n")

    cat("\n")
    print(x$summary_chain)

    cat("---\n")
    print(x$summary_exp)

    cat("---\n")
    print(x$waic)
    print(x$loo)

    invisible()
  }
