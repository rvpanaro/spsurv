rm(list = ls(all = T))
options(scipen = 9999)

require(spsurv)
citation('spsurv')
data('larynx')

Surv(larynx$time,larynx$delta)
library(xtable)
xtable(larynx)

head(larynx)

larynx$stage <-  as.factor(larynx$stage)

### Cox's Proportional Hazards Model from survival package

fitcox<- survival::coxph(Surv(time, delta) ~ age + stage, data = larynx)
fitcox

fitcox$coefficients

# 1. Proportional odds Bayesian estimates
fit1be <-spbp(Surv(time, delta) ~ age + stage, data = larynx,
              approach = "bayes", model = "po", chains = 1)
#print(fit1be, pars = "beta")
fitoddsbe <- apply(rstan::extract(fit1be, pars = "beta")$beta, 2, mean)

cbind(fitcox$coefficients, fitoddsbe)
cat( "Absolute Percentage Error%\n",
     abs(fitoddsbe-fitcox$coefficients)/fitcox$coefficients *100)

# 2. Proportional hazards Bayesian estimates
fit2be <-spbp(Surv(time, delta) ~ age + stage, data = larynx,
              approach = "bayes", model = "ph", chains = 1)
#print(fit2be, pars = "beta")
fithazardsbe <- apply(rstan::extract(fit2be, pars = "beta")$beta, 2, mean)

cbind(fitcox$coefficients, fithazardsbe)
cat( "Absolute Percentage Error%\n",
     abs(fithazardsbe-fitcox$coefficients)/fitcox$coefficients *100)

# 3. Proportional odds maximum likehlihood estimates
fit1mle <- spbp(Surv(time, delta) ~ age + stage, data = larynx,
                approach = "mle", model = "po")
#print(fit1mle)
fitoddsmle <- fit1mle$par[1:4]

cbind(fitcox$coefficients, fitoddsmle)
cat( "Absolute Percentage Error%\n",
     abs(fitoddsmle-fitcox$coefficients)/fitcox$coefficients *100)

# 4. Proportional hazards maximum likelihood estimates
fit2mle <- spbp(Surv(time, delta) ~ age + stage, data = larynx,
                approach = "mle", model = "ph")
#print(fit2mle)
fithazardsmle <- fit2mle$par[1:4]

cbind(fitcox$coefficients, fithazardsmle)
cat( "Absolute Percentage Error%\n",
     abs(fithazardsmle-fitcox$coefficients)/fitcox$coefficients *100)



