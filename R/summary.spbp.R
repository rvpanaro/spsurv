#' Bernstein Polynomial based Regression Object
#'
#' @export
#' @method summary spbp

summary.spbp <- function(spbp, interval = 0.95, ...){
    ## mle approach
    if(spbp$call$approach == "mle"){
        beta <- spbp$coefficients[1:spbp$q]
        var <- spbp$var[1:spbp$q, 1:spbp$q]

        ### Error handling ###
        # Null model
        if (is.null(spbp$coefficients)) return(spbp)

        beta2 <- beta[!(is.na(beta))] #non-missing coefs
        if(is.null(beta) | is.null(var )) stop("Input is not valid")

        se <- suppressWarnings(sqrt(diag(spbp$var)[1:spbp$q]))

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
        if (interval) {
            z <- qnorm((1 + interval)/2, 0, 1)
            output$interval <- cbind(exp(beta), exp(-beta), exp(beta - z * se),
                         exp(beta + z * se))
            dimnames(output$interval) <- list(names(beta), c("exp(coef)", "exp(-coef)",
                                                 paste("lower .", round(100 * interval, 2), sep = ""),
                                                 paste("upper .", round(100 * interval, 2), sep = "")))
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


        class(output) <- switch (spbp$call$model, "po"  = "summary.bppo.mle",
                                             "ph"  = "summary.bpph.mle",
                                             "aft" = "summary.bpaft.mle")
    }
    else{
        status <- eval(spbp$call$data)$status
        output <- list(nevent = sum(status),
                       n = length(status),
                       call = spbp$call)

        output$coef_names <- colnames(model.matrix(spbp))
        output$summarise <- rstan::summary(spbp$stanfit,
                                  probs = c((1-interval)/2, .5, interval + (1-interval)/2), pars = "beta")$summary
        exp_samp <- coda::mcmc(exp(rstan::extract(spbp$stanfit, "beta")$beta))

        output$Coefmat <- cbind(spbp$pmode$par[1:length( output$coef_names)],
                                coda::HPDinterval(log(exp_samp), prob = interval),
                                output$summarise)

        colnames(output$Coefmat) <- c("mode", "lowerHPD", "upperHPD",
                                      colnames(output$summarise))
        rownames(output$Coefmat) <- output$coef_names
        #####

        output$Coefmat2 <- cbind(exp(output$Coefmat[,1]),
                                 exp(output$Coefmat[,2]),
                                 apply(exp_samp, 2, mean),
                                 apply(exp_samp, 2, median),
                                 apply(exp_samp, 2, sd),
                                 coda::HPDinterval(exp_samp, prob = interval)
                                 )
        rownames(output$Coefmat2) <- rownames(output$Coefmat)
        colnames(output$Coefmat2) <-  c("exp(mode)", "exp(mean)", "mean_exp", "median_exp", "sd_exp",
                                        "lowerHPD_exp", "upperHPD_exp")

        output$waic <- spbp$waic
        output$loo <- spbp$loo

        class(output) <- switch (spbp$call$model, "po"  = "summary.bppo.bayes",
                                 "ph"  = "summary.bpph.bayes",
                                 "aft" = "summary.bpaft.bayes")

    }
    return(output)
}
