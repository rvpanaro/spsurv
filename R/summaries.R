summary.spbp <- function(spbp, conf.int = 0.95, scale = 1, ...) {

    if(spbp$approach == 0){
        beta <- spbp$coefficients * scale

        # Null model
        if (is.null(spbp$coefficients)) {   # Null model
            return(object)
        }
        nabeta <- !(is.na(beta))          #non-missing coefs
        beta2 <- beta[nabeta]
        if(is.null(beta) | is.null(spbp$var))
            stop("Input is not valid")
        se <- sqrt(diag(spbp$var)) * scale
        if (!is.null(spbp$naive.var)) nse <- sqrt(diag(spbp$naive.var))

        output <- list(call=spbp$call,
                       fail=spbp$fail,
                       na.action=spbp$na.action,
                       n=spbp$n,
                       loglik=spbp$loglik)

        if (!is.null(spbp$nevent))
            output$nevent <- spbp$nevent

        if (is.null(spbp$naive.var)) {
            tmp <- cbind(beta, exp(beta), se, beta/se,
                         pchisq((beta/ se)^2, 1, lower.tail=FALSE))
            dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)",
                                                 "se(coef)", "z", "Pr(>|z|)"))
        }
        else {
            tmp <- cbind(beta, exp(beta), nse, se, beta/se,
                         pchisq((beta/ se)^2, 1, lower.tail=FALSE))
            dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)",
                                                 "se(coef)", "robust se", "z", "Pr(>|z|)"))
        }
        output$coefficients <- tmp

        if (conf.int) {
            z <- qnorm((1 + conf.int)/2, 0, 1)
            tmp <- cbind(exp(beta), exp(-beta), exp(beta - z * se),
                         exp(beta + z * se))
            dimnames(tmp) <- list(names(beta), c("exp(coef)", "exp(-coef)",
                                                 paste("lower .", round(100 * conf.int, 2), sep = ""),
                                                 paste("upper .", round(100 * conf.int, 2), sep = "")))
            output$conf.int <- tmp
        }

        df <- length(beta2)
        logtest <- -2 * (spbp$loglik[1] - spbp$loglik[2])
        output$logtest <- c(test=logtest,
                            df=df,
                            pvalue= pchisq(logtest, df, lower.tail=FALSE))
        output$sctest <- c(test=spbp$score,
                           df=df,
                           pvalue= pchisq(spbp$score, df, lower.tail=FALSE))
        output$rsq<-c(rsq=1-exp(-logtest/spbp$n),
                      maxrsq=1-exp(2*spbp$loglik[1]/spbp$n))
        output$waldtest<-c(test=as.vector(round(spbp$wald.test, 2)),
                           df=df,
                           pvalue= pchisq(as.vector(spbp$wald.test), df,
                                          lower.tail=FALSE))
        if (!is.null(spbp$rscore))
            output$robscore<-c(test=spbp$rscore,
                               df=df,
                               pvalue= pchisq(spbp$rscore, df, lower.tail=FALSE))
        output$used.robust<-!is.null(spbp$naive.var)

        if (!is.null(spbp$concordance)) {
            # throw away the extra info, in the name of backwards compatability
            output$concordance <- spbp$concordance[6:7]
            names(output$concordance) <- c("C", "se(C)")
        }
        ifelse(spbp$model == 1,
               class(output) <- "summary.bpph.mle",
               class(output) <- "summary.bppo.mle")
    }
    else{
        samp <- coda::mcmc(rstan::extract(spbp$stanfit, pars = "beta")$beta)
        output <- list(summary = cbind(Median = apply(samp, 2, median), summary(samp)[[1]],
                                       HPDL = coda::HPDinterval(samp)[,1], HPDU = coda::HPDinterval(samp)[,2])[,c(1,2,3,6,7)])
        # rownames(output$summary) <- colnames(Z)

        ifelse(spbp$model == 1,
               class(output) <- "summary.bpph.bayes",
               class(output) <- "summary.bppo.bayes")
    }
    return(output)
}
