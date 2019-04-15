rm(list = ls(all = T))
options(scipen = 9999)

require(spsurv)
citation('spsurv')
data('larynx')

head(larynx)

larynx$stage <-  as.factor(larynx$stage)
# spsurv::bpph(Surv(time, delta) ~ age + stage, data = larynx)

fit <-spbp(Surv(time, delta) ~ age + stage, data = larynx, approach = "bayes")
# spsurv::bpph(Surv(time, delta) ~ age + stage, data = larynx)

fit <- spbp(Surv(time, delta) ~ age + stage, data = larynx, approach = "mle")
fit

rstan::traceplot(fit)

