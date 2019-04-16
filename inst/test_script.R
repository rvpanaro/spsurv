rm(list = ls(all = T))
options(scipen = 9999)

require(spsurv)
citation('spsurv')
data('larynx')

head(larynx)

larynx$stage <-  as.factor(larynx$stage)

fit1be <-spbp(Surv(time, delta) ~ age + stage, data = larynx,
              approach = "bayes", baseline = "po", chains = 1)
print(fit1be)

fit2be <-spbp(Surv(time, delta) ~ age + stage, data = larynx,
              approach = "bayes", baseline = "ph", chains = 1)
print(fit2be)

fit1mle <- spbp(Surv(time, delta) ~ age + stage, data = larynx,
                approach = "mle", model = "po", init = 0)
print(fit1mle)

fit2mle <- spbp(Surv(time, delta) ~ age + stage, data = larynx,
                approach = "mle", model = "ph", init = 0)
print(fit2mle)




