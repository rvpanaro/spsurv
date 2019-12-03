#' Bernstein Polynomial based Regression Object
#'
#' @export
#' @method print spbp

print.spbp <-
  function(spbp, digits = max(getOption('digits')-4, 3),
           signif.stars = getOption("show.signif.stars"), ...){

  savedig <- options(digits = digits)
  on.exit(options(savedig))

  if (!is.null(spbp$call)) {
    cat("\n---")
    cat("\n")
    cat("Call:\n")
    dput(spbp$call)
    cat("\n")
  }

  if(spbp$call$approach == "mle"){

    coef <- as.array(spbp$coefficients[1:spbp$q])
    var <- as.array(spbp$var[1:spbp$q, 1:spbp$q])

    ### Error handling ###
    # Null model
    if (is.null(spbp$coefficients)) return(spbp)

    coef2 <- coef[!(is.na(coef))] #non-missing coefs
    if(is.null(coef) | is.null(var )) stop("Input is not valid")

    se <- as.array(suppressWarnings(sqrt(diag(spbp$var)[1:spbp$q])))

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
    summarise <- rstan::summary(spbp$stanfit, pars = "beta_std")$summary
    design <- as.matrix(model.matrix(spbp))
    p <- ncol(design)
    Coef <- cbind(matrix(summarise[, 1], nrow = p),
                  matrix(coda::HPDinterval(coda::mcmc(rstan::extract(spbp$stanfit, "beta_std")$beta_std)), nrow = p),
                  matrix(summarise[, -c(1, 5, 7, 9, 10)], nrow = p))

    rownames(Coef) <-  colnames(design)
    colnames(Coef) <- c("mean", "lowerHPD", "upperHPD", colnames(summarise)[-c(1, 5, 7, 9, 10)])
    print(Coef, digits = digits)

    cat("---\n")
    cat("\n WAIC Estimate= ", sprintf('%.3f', spbp$waic$estimates[3,1]))
    cat("      WAIC SE= ", sprintf('%.3f', spbp$waic$estimates[3,2], "\n"))
    cat("\n LOOIC Estimate= ", sprintf('%.3f', spbp$loo$estimates[3,1]))
    cat("      LOOIC SE= ", sprintf('%.3f', spbp$loo$estimates[3,2], "\n"))
    cat("\n---\n")
    cat("\n")
  }
}

