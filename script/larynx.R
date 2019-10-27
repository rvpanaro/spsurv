rm(list = ls())
library("KMsurv")
data("larynx")

library("spsurv")
fitmle <- spbp(Surv(time, delta)~ age + factor(stage),
     data = larynx)

fitmle

fitbe <- spbp(Surv(time, delta)~ age + factor(stage), scale = T,
               data = larynx, approach  = "bayes")

fitbe

traceplot(fitbe$stanfit, pars = c("beta_std", "gamma_std"))
stan_dens(fitbe$stanfit, pars = c("gamma_std", "gamma_std"))

smle1 <- survivor(fitmle, newdata = data.frame(77, 0, 0, 0))
sbe1 <- survivor(fitbe, newdata = data.frame(77, 0, 0, 0))
cbind(smle1, sbe1)

smle1 <- survivor(fitmle, newdata = data.frame(77, 1, 0, 0))
sbe1 <- survivor(fitbe, newdata = data.frame(77, 1, 0, 0))

smle2 <- survivor(fitmle, newdata = data.frame(77, 0, 1, 0))
sbe2 <- survivor(fitbe, newdata = data.frame(77, 0, 1, 0))

smle3 <- survivor(fitmle, newdata = data.frame(77, 0, 0, 1))
sbe3 <- survivor(fitbe, newdata = data.frame(77, 0, 0, 1))

plot(smle3, col = 2)
points(sbe3, col = 4)

