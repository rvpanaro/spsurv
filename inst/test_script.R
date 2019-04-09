rm(list = ls(all = T))

require(spsurv)
citation('spsurv')
data('larynx')

head(larynx)

larynx$stage <-  as.factor(larynx$stage)
# spsurv::bpph(Surv(time, delta) ~ age + stage, data = larynx)

fit <- bpph(Surv(time, delta) ~ age + stage, data = larynx, approach = "mle")
names(fit)


fit <- spsurv::bpph(Surv(time, delta) ~ age + stage, data = larynx, approach = "bayes")
# spsurv::bpph(Surv(time, delta) ~ age + stage, data = larynx)

rstan::traceplot(fit)


prior <- list(a_gamma = 0, b_gamma = 0, m_beta = 0, S_beta = 0)
standata <- do.call(c, list(prior, prior))
