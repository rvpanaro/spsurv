source(file = 'inst/guidelines_load.R')
# library(spsurv)
data("larynx"); str(larynx)

fit1 <- spbp(Surv(time, delta) ~ age + factor(stage),
            approach = 'mle', model = 'ph', data = larynx)
summary(fit1)
fit1$stanfit$return_code

fit2 <- spbp(Surv(time, delta) ~ age + factor(stage),
     approach = 'bayes', model = 'aft', data = larynx, chain = 1)

class(fit2)
summary(fit2)
fit2$stanfit$return_code
getwd()
