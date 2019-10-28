getwd()
setwd("../MC")
# source('script/load.R')

rm(list = ls())
try(dev.off())
?spsurv

library(spsurv)
library("data.table")
dat <- fread("../MC/aft_llogis/data_aft_llogis_899.txt")
head(dat)
mean(dat$status)

## BPPH MLE
fitmle <- spbp(Surv(y, status)~ ., data = dat,
               approach = "mle", model = "aft")
## BPPH BE
fitbe <- spbp(Surv(y, status)~., data = dat,
               approach = "bayes", model = "aft")
fitbe
new <- data.frame(numerical = c(mean(dat$numerical)),
                  categorical = c(0,1))

plot(sort(fitmle$y[,1]), survivor.spbp(fitmle, newdata = new))
plot(sort(fitmle$y[,1]), survivor.spbp(fitmle, newdata = new), pch = 19, ylim = c(0,1), bty = 'n')
points(sort(fitbe$y[,1]), apply(survivor.spbp(fitbe, newdata = new), 2, mean), pch = 4)
## actual survival
eta <- as.matrix(new) %*% c(2,-1)
scaleT <- 8
shape <- 2
location = log(scaleT) - shape * eta ## logis location parameter
scale <- exp(-location/shape)

library("eha")
curve(pllogis(x, shape = shape, scale = scale[1,1], lower = F), add = T, lwd = 2, lty = 2, col = "darkgrey")
curve(pllogis(x, shape = shape, scale = scale[2,1], lower = F), add = T, lwd = 2, lty = 2, col = "darkgrey")
legend("bottomleft", c("Actual (dashed)", "Bayes", "MLE", "Breslow"),
       col = c("darkgrey", "black", "black", "black"),
       pch = c(95, 4, 19,95), bty = "n")
