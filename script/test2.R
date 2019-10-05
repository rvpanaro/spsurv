source(file = "inst/load.R")
library(spsurv)

# weibull aft
db1 <- sim_surv(100, model = "aft")
fit1 <- spbp(Surv(y,status) ~ x1 + x2,  model = "aft",
             approach = "mle", data = db1)
t1 <- Sys.time()
fit1_bayes <- spbp(Surv(y,status) ~ x1 + x2,  model = "aft",
             approach = "bayes", data = db1)
t2 <- Sys.time()
t2 -t1
fit1_bayes

#
# pred1 <- fit1$linear.predictors

## weibull ph

db2 <- sim_surv(100)
fit2 <- spbp(Surv(y,status)~ x1 + x2,  model = "ph",
             approach = "mle", data = db2)
fit2
summary(fit2)

fit2_bayes <- spbp(Surv(y,status)~ x1 + x2,  model = "ph",
             approach = "bayes", data = db2, chains = 1)
fit2_bayes
summary(fit2_bayes)

library(loo)
extract_log_lik(fit2_bayes$stanfit)

surv_ref2 <- function(x)pweibull(x, shape = gamma,
                                       scale = (lambdaT*exp(new[-1] %*% c(beta1, beta2))^(-1/gamma)), lower = FALSE)

## weibull po
db3 <- sim_surv(100, model = "po", scaleC = .4)
fit3 <- spbp(Surv(y,status)~ x1 + x2,  model = "po",
             approach = "mle", data = db3)
pred3 <- fit3$linear.predictors
surv_ref3 <- function(x)pweibull(x, shape = gamma,
                                       scale = (lambdaT*exp(new[-1] %*% c(beta1, beta2))^(-1/gamma)), lower = FALSE)

## llogis aft
db4 <- sim_surv(100, dist = "llogis", model = "aft")
fit4 <- spbp(Surv(y,status)~ x1 + x2,  model = "aft",
             approach = "mle", data = db4)
pred4 <- fit4$linear.predictors
surv_ref4 <- function(x)pweibull(x, shape = gamma,
                                       scale = (lambdaT*exp(new[-1] %*% c(beta1, beta2))^(-1/gamma)), lower = FALSE)

## llogis ph
db5 <- sim_surv(100, dist = "llogis", model = "ph", scaleT = .3)
fit5 <- spbp(Surv(y,status)~ x1 + x2,  model = "ph",
             approach = "mle", data = db5)
pred5 <- fit5$linear.predictors
surv_ref5 <- function(x)pweibull(x, shape = gamma,
                                       scale = (lambdaT*exp(new[-1] %*% c(beta1, beta2))^(-1/gamma)), lower = FALSE)

## llogis po
db6 <- sim_surv(100, dist = "llogis", model = "po")
fit6 <- spbp(Surv(y,status)~ x1 + x2,  model = "po",
             approach = "mle", data = db6)
pred6 <- fit6$linear.predictors
surv_ref6 <- function(x)pweibull(x, shape = gamma,
                                       scale = (lambdaT*exp(new[-1] %*% c(beta1, beta2))^(-1/gamma)), lower = FALSE)

