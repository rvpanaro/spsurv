source(file = 'inst/guidelines_load.R')

n = 100
beta1 = 2; beta2 = -1
lambdaT = .002 # baseline hazard
lambdaC = .004  # hazard of censoring

x1 = rnorm(n,0)
x2 = rnorm(n,0)
# true event time
t = rweibull(n, shape=1, scale = lambdaT*exp(-beta1*x1-beta2*x2))
c = rweibull(n, shape=1, scale = lambdaC)   #censoring time
time = pmin(t,c)  #observed time is min of censored and true
event = time == t   # set to 1 if event is observed

X <- data.frame(time, as.numeric(event), x1, x2)

library(survival) #standard model
fit_survival <- survreg(Surv(time, event)~ x1 + x2)
summary(fit_survival)


library(spsurv) #bernstein polynomial based regression
fit_spbp <- spbp(Surv(time, event) ~ x1 + x2, data = X,
            model = 'ph', approach = 'mle')

plot(fit_spbp)
