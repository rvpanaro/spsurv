## spbp object printings
print.spbp <- function(spbp, digits = max(getOption('digits')-3, 3),
                       signif.stars = getOption("show.signif.stars"), ...){

  savedig <- options(digits = digits)
  on.exit(options(savedig))

  if(spbp$approach == "mle"){

    coef <- spbp$coefficients
    var <- spbp$var

    ### Error handling ###
    # Null model
    if (is.null(spbp$coefficients)) return(spbp)

    coef2 <- coef[!(is.na(coef))] #non-missing coefs
    if(is.null(coef) | is.null(var )) stop("Input is not valid")

    se <- sqrt(diag(spbp$var))

    Coefmat  <- cbind(coef, exp(coef), se, coef/se,
                                  pchisq((coef/ se)^2, 1, lower.tail=FALSE))
    dimnames(Coefmat) <- list(names(coef), c("coef", "exp(coef)",
                                                         "se(coef)", "z", "Pr(>|z|)"))


    if (!is.null(spbp$call)) {
      cat("Call:\n")
      dput(spbp$call)
      cat("\n")
    }

    if(!is.null(spbp$coefficients)) {
      cat("\n")
      printCoefmat(Coefmat, digits = digits,
                   signif.stars = signif.stars, ...)
    }
    if(!is.null(spbp$loglik)) {
      cat("\n Loglik(model)= ", spbp$loglik[2])
      cat("      Loglik(no predictors)= ", spbp$loglik[1], "\n")
    }
    logtest <- -2 * (spbp$loglik[1] - spbp$loglik[2])
    cat("      Chisq= ", logtest," on ", spbp$q, " degrees of freedom ",
        pchisq(logtest, spbp$q, lower.tail=FALSE), "\n")

    if(!is.null(spbp$n)) {
        cat("n= ", spbp$n)
    }
  }
  else{
    print('no print methos yet')
  }
}

## summary printings
print.summary.spbp.mle <- function(x, digits = max(getOption('digits')-3, 3),
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

print.summary.bppo.mle <- function(...){
  cat("Bernstein Polynomial based Proportional Odds model\n")
  print.summary.spbp.mle(...)
}

print.summary.bpph.mle <- function(...){
  cat("Bernstein Polynomial based Proportional Hazards model\n")
  print.summary.spbp.mle(...)
}

print.summary.bpaft.mle <- function(...){
  cat("Bernstein Polynomial based Accelerated Failure Time model\n")
  print.summary.spbp.mle(...)
}

