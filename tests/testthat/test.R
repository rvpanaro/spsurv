
library(spsurv)
data("veteran")
library("KMsurv")
data("larynx")

# ML approach:
fit <- spbp(Surv(time, delta)~age+factor(stage),
            approach = "mle",  data = larynx)
summary(fit)
vcov(fit)
coef(fit)
model.matrix(fit)
# stan_dens(fit)
# traceplot(fit)
# extract(fit)

# Bayesian approach:
fit2 <- spbp(Surv(time, delta) ~ age + factor(stage),
             approach = "bayes",  data = larynx,
             iter = 20, chains = 1, warmup = 10)
summary(fit2)
# vcov(fit2)
# coef(fit2)
model.matrix(fit2)
stan_dens(fit2)
traceplot(fit2)
extract(fit2)

## PH model
bpph(Surv(time, delta)~age+factor(stage),  data = larynx)## PO model
bpph(Surv(time, delta)~age+factor(stage), approach = "bayes",  data = larynx,
     iter = 20, chains = 1, warmup = 10)## PH model

## PO model
bppo(Surv(time, delta)~age+factor(stage),  data = larynx)## PO model
bppo(Surv(time, delta)~age+factor(stage),  approach = "bayes",  data = larynx,
     iter = 20, chains = 1, warmup = 10)## PO model

## AFT model
bpaft(Surv(time, delta)~age+factor(stage),  data = larynx)## AFT model
bpaft(Surv(time, delta)~age+factor(stage),  approach = "bayes",  data = larynx,
     iter = 20, chains = 1, warmup = 10)## PO model
