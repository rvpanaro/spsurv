summary.spbp <- function(spbp, conf.int = 0.95, ...) {
    ## mle approach

    if(spbp$approach == "mle"){
        beta <- spbp$coefficients[1:spbp$q]
        var <- spbp$var[1:spbp$q, 1:spbp$q]

        ### Error handling ###
        # Null model
        if (is.null(spbp$coefficients)) return(spbp)

        beta2 <- beta[!(is.na(beta))] #non-missing coefs
        if(is.null(beta) | is.null(var )) stop("Input is not valid")

        se <- sqrt(diag(spbp$var)[1:spbp$q])

        output <- list(call = spbp$call,
                       return_code = spbp$return_code,
                       n = spbp$n,
                       loglik = spbp$loglik)

        if (!is.null(spbp$nevent))
         output$nevent <- spbp$nevent

        output$coefficients  <- cbind(beta, exp(beta), se, beta/se,
                         pchisq((beta/ se)^2, 1, lower.tail=FALSE))
        dimnames(output$coefficients) <- list(names(beta), c("coef", "exp(coef)",
                                                 "se(coef)", "z", "Pr(>|z|)"))
        if (conf.int) {
            z <- qnorm((1 + conf.int)/2, 0, 1)
            output$conf.int <- cbind(exp(beta), exp(-beta), exp(beta - z * se),
                         exp(beta + z * se))
            dimnames(output$conf.int) <- list(names(beta), c("exp(coef)", "exp(-coef)",
                                                 paste("lower .", round(100 * conf.int, 2), sep = ""),
                                                 paste("upper .", round(100 * conf.int, 2), sep = "")))
        }
        df <- length(beta2)
        logtest <- -2 * (spbp$loglik[1] - spbp$loglik[2])
        output$logtest <- c(test=logtest,
                            df=df,
                            pvalue= pchisq(logtest, df, lower.tail=FALSE))
        output$rsq<-c(rsq=1-exp(-logtest/spbp$n),
                      maxrsq=1-exp(2*spbp$loglik[1]/spbp$n))
        output$waldtest<-c(test=as.vector(round(spbp$wald.test, 2)),
                           df=df,
                           pvalue= pchisq(as.vector(spbp$wald.test), df,
                                          lower.tail=FALSE))


        class(output) <- switch (spbp$model, "po"  = "summary.bppo.mle",
                                             "ph"  = "summary.bpph.mle",
                                             "aft" = "summary.bpaft.mle")
    }
    else{
        samp <- coda::mcmc(rstan::extract(spbp$stanfit, pars = "beta")$beta)
        output <- list(summary = cbind(Median = apply(samp, 2, median), summary(samp)[[1]],
                                       HPDL = coda::HPDinterval(samp)[,1], HPDU = coda::HPDinterval(samp)[,2])[,c(1,2,3,6,7)])
        # rownames(output$summary) <- colnames(Z)
        class(output) <- switch (spbp$model, "po"  = "summary.bppo.bayes",
                                             "ph"  = "summary.bpph.bayes",
                                             "aft" = "summary.bpaft.bayes")
    }
    return(output)
}
