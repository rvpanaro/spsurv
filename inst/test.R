source('inst/load.R')
rm(list = ls())
## Weibull simulated data: Inverse transform method

n <- 100
beta1 <- 2; beta2 = -1
lambdaT <- .002 # baseline hazard
lambdaC <- .004  # hazard of censoring

x1 <- rnorm(n, 0)
x2 <- rnorm(n, 0)

gamma <- 2

t <- rweibull(n, shape = gamma, scale = (lambdaT*exp(beta1*x1 + beta2*x2)^(-1/gamma))) #proportional hazards
c <- rweibull(n, shape = gamma, scale = lambdaC)   #censoring time
time <- pmin(t,c)  #observed time is min between censored and event
status <- as.numeric(time == t)   # set to 1 if event is observed
mean(status)

dat <- data.frame(time, status, x1, x2)

library(spsurv) # semiparametric survival

fit <- spbp(Surv(time, status) ~ x1 + x2 + frailty.gamma(x1), data = dat,
            model = 'ph', approach = 'bayes', chain = 1)
View(fit)
print(fit$stanfit)
summary(fit)

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
