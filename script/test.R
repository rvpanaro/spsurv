source('inst/load.R')
rm(list = ls())
## Weibull simulated data: Inverse transform method


library(spsurv) # semiparametric survival
data("llogis_aft")

?sim_weibull()
fitph <- spbp(Surv(time, status) ~ x1 + x2 + frailty.gamma(x1), data = dat,
            model = 'ph', approach = 'bayes', chain = 1)

summary(fitph)

cox <- coxph(Surv(time, status) ~ x1 + x2, data = dat)
cox
summary(cox)

round(diag(fit$var), 3)
data("lung")

# Random institutional effect
cox_fra <- coxph(formula = Surv(time, status) ~ I(age, 6) + frailty(inst, df=4), lung)
model.matrix(cox_fra)
View(formula)
survival:::coxph
