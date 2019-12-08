rm(list = ls())

library("spsurv")

## 1. Larynx
library("KMsurv")
data("larynx")

# BPPH
bpphmle1 <- spbp(formula = Surv(time, delta) ~ age + factor(stage),
     data = larynx,
     approach = "mle",
     model = "ph")

summary(bpphmle1)

bpphbe1 <- spbp(Surv(time, delta)~ age + factor(stage),
               data = larynx,
               approach  = "bayes",
               model = "ph")
summary(bpphbe1)

# BPPO
bppomle1 <- spsurv::spbp(formula = Surv(time, delta)~ age + factor(stage),
                         data = larynx,
                         approach = "mle",
                         model = "po")
summary(bppomle1)

bppobe1 <- spbp(Surv(time, delta)~ age + factor(stage), scale = T,
                data = larynx,
                approach  = "bayes",
                model = "po")
summary(bppobe1)

# BPAFT
bpaftmle1 <- spsurv::spbp(formula = Surv(time, delta)~ age + factor(stage),
                         data = larynx,
                         approach = "mle",
                         model = "aft")
summary(bpaftmle1)

bpaftbe1 <- spbp(Surv(time, delta)~ age + factor(stage), scale = T,
                data = larynx,
                approach  = "bayes",
                model = "aft")
summary(bpaftbe1)

## 2. Lung
data("veteran")
attach(veteran)

veteran$celltype <- factor(celltype, c("large", "adeno", "smallcell", "squamous"))
fitmle2 <- spsurv::spbp(formula = Surv(time, status) ~ karno + celltype,
                       data = veteran, approach = "mle", model = "aft", scale = T)
summary(fitmle2)

fitbe2 <- spsurv::spbp(formula = Surv(time, status) ~ karno + celltype,
                        data = veteran, approach = "bayes", model = "po")
summary(fitbe2)

library("timereg")
prop.odds(formula = Event(time, status) ~ karno + celltype,
           data = veteran)

library("survival")
svg <- survreg(formula = Surv(time, status) ~ karno + celltype,
          data = veteran, dist = "loglogistic")
summary(svg)

confint(svg)
