# library(spsurv)
library(survival)
data("veteran")
str(veteran)

fit <- spbp(Surv(time, status) ~ karno + celltype,
     approach = 'mle', model = 'aft', data = veteran)

names(fit)
fit$stanfit
fit$model
fit$approach
summary(fit)
class(fit)
traceplot.spbp(fit, pars = 'beta')
traceplot(fit, pars = 'beta')
survfit(fit)

