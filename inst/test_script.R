rm(list = ls(all = T))
require(spsurv)
citation('spsurv')
data('larynx')

head(larynx)

larynx$stage <-  as.factor(larynx$stage)
# spsurv::bpph(Surv(time, delta) ~ age + stage, data = larynx)

fit <- spsurv::bpph(Surv(time, delta) ~ age + stage, data = larynx, approach = "frequentist")
# spsurv::bpph(Surv(time, delta) ~ age + stage, data = larynx)

