# source(file = 'inst/guidelines_load.R')

rm(list = ls())
options(scipen = 9999)

library(spsurv)
library(survival)
data("larynx")
str(larynx)

fit <- spbp(Surv(time, delta) ~  factor(stage) + age,
             approach = 'mle', model = 'ph', data = larynx)

fit
plot(fit)
summary(fit)
fit$call$data
plot(fit)

fit_survival <- coxph(Surv(time, delta) ~ factor(stage) + age, data = larynx)
fit_survival$coefficients
plot(survfit(fit_survival))
