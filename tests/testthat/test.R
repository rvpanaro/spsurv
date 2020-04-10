library(spsurv)
data("veteran")

#---- bpfits.R ----

## PH model
ph_fit <- bpph(Surv(time, status)~karno+factor(celltype),  data = veteran)## PO model
print(ph_fit)
survivor(ph_fit)
summary(ph_fit)
residuals(ph_fit)
print(summary(ph_fit))
expect_warning(traceplot(ph_fit))
expect_warning(extract(ph_fit))
expect_warning(stan_dens(ph_fit))

ph_fit <- bpph(Surv(time, status)~karno+factor(celltype), approach = "bayes",  data = veteran,
     iter = 10, chains = 1, cores = 1)## PH model
print(ph_fit)
survivor(ph_fit)
summary(ph_fit)
print(summary(ph_fit))

expect_warning(confint(ph_fit))
expect_warning(coef(ph_fit))
expect_warning(vcov(ph_fit))

## PO model
po_fit <- bppo(Surv(time, status)~karno+factor(celltype),  data = veteran)## PO model
print(po_fit)
summary(po_fit)
survivor(po_fit)
print(summary(po_fit))

po_fit <- bppo(Surv(time, status)~karno+factor(celltype),  approach = "bayes",  data = veteran,
     iter = 10, chains = 1, cores = 1)## PO model
print(po_fit)
summary(po_fit)
survivor(po_fit)
print(summary(po_fit))

## AFT model
aft_fit <- bpaft(Surv(time, status)~karno+factor(celltype),  data = veteran)## AFT model
print(aft_fit)
summary(aft_fit)
survivor(aft_fit)
print(summary(aft_fit))

aft_fit <- bpaft(Surv(time, status)~karno+factor(celltype),  approach = "bayes",  data = veteran,
      iter = 10, chains = 1, cores = 1)## PO model
print(aft_fit)
summary(aft_fit)
survivor(aft_fit)
print(summary(aft_fit))

#--- rsurv ---
rsurv(100)
rsurv(100, dist = "log-logistic")

#---- handlers.R ----

expect_error(spbp(formula = Surv(time, status) ~ karno + factor(celltype),
     degree = 1/2,
     data = veteran))

expect_error(spbp(formula = Surv(time, status) ~ karno + factor(celltype),
     data = veteran, chaind = 1))

expect_error(spbp(formula = time ~ karno + factor(celltype),
     data = veteran))

time2 <- veteran$time
expect_error(spbp(formula = Surv(time = time, time2 = time2, status) ~ karno + factor(celltype),
     data = veteran))

expect_error(spbp(formula = Surv(type = "left", time = time, status) ~ karno + factor(celltype),
     data = veteran))

expect_error(spbp(formula = Surv(time = time, status) ~ status+ karno,
     data = veteran))

#---- spbp.R

# ML approach:
fit <- spbp(Surv(time, status)~karno+factor(celltype),
            approach = "mle",  data = veteran)
summary(fit)
vcov(fit)
coef(fit)
model.matrix(fit)
print(fit)
# stan_dens(fit)
# traceplot(fit)
# extract(fit)

# Bayesian approach:
fit2 <- spbp(Surv(time, status) ~ karno + factor(celltype),
             approach = "bayes",  data = veteran, chains = 1, iter = 10, cores = 1)
print(fit2)
summary(fit2)
# vcov(fit2)
# coef(fit2)
model.matrix(fit2)
stan_dens(fit2)
traceplot(fit2)
extract(fit2)


