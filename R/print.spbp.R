print.spbp <-
  function(spbp, digits = max(getOption('digits')-3, 3),
           signif.stars = getOption("show.signif.stars"), ...){

  savedig <- options(digits = digits)
  on.exit(options(savedig))

  if (!is.null(spbp$call)) {
    cat("Call:\n")
    dput(spbp$call)
    cat("\n")
  }

  if(spbp$call$approach == "mle"){

    coef <- spbp$coefficients[1:spbp$q]
    var <- spbp$var[1:spbp$q, 1:spbp$q]

    ### Error handling ###
    # Null model
    if (is.null(spbp$coefficients)) return(spbp)

    coef2 <- coef[!(is.na(coef))] #non-missing coefs
    if(is.null(coef) | is.null(var )) stop("Input is not valid")

    se <- suppressWarnings(sqrt(diag(spbp$var[1:spbp$q, 1:spbp$q])))

    Coefmat  <- cbind(coef, exp(coef), se, coef/se,
                                  pchisq((coef/ se)^2, 1, lower.tail=FALSE))
    dimnames(Coefmat) <- list(names(coef), c("coef", "exp(coef)",
                                                         "se(coef)", "z", "Pr(>|z|)"))

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
    cat("\n")
    summarise <- rstan::summary(spbp$stanfit, pars = "beta")$summary
    Coef <- cbind(summarise[, 1],
                  exp(summarise[,1]),
                  coda::HPDinterval(coda::mcmc(rstan::extract(spbp$stanfit, "beta")$beta)),
                  summarise[, -c(1, 5, 7, 9, 10)])
    rownames(Coef) <-  all.vars(spbp$call$formula)[-c(1,2)]
    colnames(Coef) <- c("mean", "exp(mean)", "lowerHPD", "upperHPD", colnames(summarise[, -c(1, 5, 7, 9, 10)]))
    print(Coef)

    cat("---\n")
    print(cbind(t(spbp$waic$estimates), t(spbp$loo$estimates)))
  }
}

