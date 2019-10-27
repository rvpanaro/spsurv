getwd()
setwd("../MC")
# source('script/load.R')

rm(list = ls())
dev.off()
?spsurv

library(spsurv)
library("data.table")
dat <- fread("../MC/po_weibull/data_po_weibull_111.txt")
head(dat)
mean(dat$status)

fitmle <- spbp(Surv(y, status)~ ., data = dat,
               approach = "mle", model = "po", scale = F)

summary(fitmle)

# coxph(Surv(y, status)~., data = dat)
library("timereg")
prop.odds(Event(y, status)~., data = dat)

fitbe1 <- spbp(Surv(y, status)~., data = dat,
              approach = "bayes", model = "ph",
              priors = list(beta = c("normal(0,4)", "normal(0,4)"),
                            gamma = "lognormal(0,4)")
              )
summary(fitbe1)

X11()
traceplot(fitbe1$stanfit, pars = c("beta", "gamma"))

X11()
stan_dens(fitbe1$stanfit, pars = c("beta", "gamma"))

