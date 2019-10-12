#' Bernstein Polynomial based Regression Object
#'
#' @export
#' @method print summary.spbp.mle
#'

print.summary.spbp.mle <- ## summary printings
  function(x, digits = max(getOption('digits')-4, 3),
                                   signif.stars = getOption("show.signif.stars"), ...) {
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

  if (nrow(x$coef)==0) {   # Null model
    cat ("Null model\n")
    return()
  }
  if(!is.null(x$coefficients)) {
    cat("\n")
    printCoefmat(x$coefficients, digits = digits,
                 signif.stars = signif.stars, ...)
  }
  if(!is.null(x$conf.int)) {
    cat("\n")
    print(x$conf.int)
  }
  cat("\n")
  pdig <- max(1, getOption("digits")-4)  # default it too high IMO
  cat("Likelihood ratio test= ", format(round(x$logtest["test"], 2)), "  on ",
      x$logtest["df"], " df,", "   p=",
      format.pval(x$logtest["pvalue"], digits=pdig),
      "\n", sep = "")
  cat("Wald test            = ", format(round(x$waldtest["test"], 2)), "  on ",
      x$waldtest["df"], " df,", "   p=",
      format.pval(x$waldtest["pvalue"], digits=pdig),
      "\n", sep = "")
  invisible()
}
