rm(list = ls())
?spbp
library("data.table")
dat <- fread("../MC/ph_weibull/data_ph_weibull_128.txt")

fitmle <- spbp(Surv(y, status)~., data = dat, approach = "mle")
print(fitmle)

fitbe <- spbp(Surv(y, status)~., data = dat, approach = "bayes")
print(fitbe)

class(fitmle$coefficients)
class(c(1,colMeans(extract(fitbe$stanfit)$beta)))

cox <- coxph(Surv(y, status)~., data = dat)

plot(survfit(cox)$surv, col = 1)
points(survivor(fitmle), col = 2)
points(survivor(fitbe), col = 4)
