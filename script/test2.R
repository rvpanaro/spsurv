library(spsurv)

dt <- sim_surv(n = 200, model = "aft")

fit <- survreg(Surv(y, status)~ x1 + x2, data = dt)
fit

dt <- sim_surv(n = 200, scaleC = 42, model = "aft",
               dist = "llogis")

fit <- survreg(Surv(y, status)~ x1 + x2, data = dt)
fit
