getwd()
setwd("../MC")
# source('script/load.R')

rm(list = ls())
try(dev.off())
?spsurv

library(spsurv)
library("data.table")
dat <- fread("../MC/ph_weibull/data_ph_weibull_888.txt")
head(dat)
mean(dat$status)

## CoxPH
cox <- coxph(Surv(y, status)~ ., data = dat)

## BPPH MLE
fitmle <- spbp(Surv(y, status)~ ., data = dat,
               approach = "mle", model = "ph")
## BPPH BE
fitbe <- spbp(Surv(y, status)~., data = dat,
               approach = "bayes", model = "ph")

new <- data.frame(numerical = mean(dat$numerical),
                  categorical = c(0))

coxsurv <- survfit(cox, new)
summary(coxsurv)

survmle <- survivor(fitmle, new)

plot(coxsurv, conf.int = F, lwd = 2, bty = 'n', ylim = c(0,1))
points(sort(fitmle$y[,1]), survmle, pch = 19)
points(sort(fitbe$y[,1]), apply(survivor.spbp(fitbe, newdata = new), 2, mean), pch = 4)
## actual survival
eta <- as.matrix(new) %*% c(2,-1)
scaleT <- 8
shape <- 2
scale <- (scaleT * exp(eta) ^ (-1 / shape))

curve(pweibull(x, scale = scale[1,1], shape = shape, lower = F), add = T, lwd = 2, lty = 2, col = "darkgrey")
curve(pweibull(x, scale = scale[2,1], shape = shape, lower = F), add = T, lwd = 2, lty = 2, col = "darkgrey")
legend("bottomleft", c("Actual (dashed)", "Bayes", "MLE", "Breslow"),
       col = c("darkgrey", "black", "black", "black"),
       pch = c(95, 4, 19,95), bty = "n")
