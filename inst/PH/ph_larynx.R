rm(list = ls()); graphics.off();
library(survival); library(KMsurv)
setwd("~/Documents/spsurv/inst/PH")
set.seed(1)

################################## MAXIMUM LIKELIHOOD ESTIMATION ##################################
# Cox Model
data("larynx")
dat <- larynx;
dat$status <- dat$delta

dat$age <- scale(larynx$age, scale = F)
dat$stage2 = scale(as.numeric(larynx$stage == 2), scale = F)
dat$stage3 = scale(as.numeric(larynx$stage == 3), scale = F)
dat$stage4 = scale(as.numeric(larynx$stage == 4), scale = F)

cox <- coxph(Surv(dat$time, dat$delta) ~ age + stage2 + stage3 + stage4,
             data = dat)

for (i in 1:4){
  print((cox$coefficients[i] + c(-1.96, + 1.96)*sqrt(cox$var[i, i])), digits = 4)
}

source('fun.R')

# Classical
fit1 <- optim(c(rep(1, m), rep(0, q)), fn = loglik, method = "BFGS", control = list(fnscale=-1), m = m,
              delta = delta, Z = Z, b = b, B = B, hessian = TRUE)

out1 <- list(coefficients=c(age = fit1$par[11], stage2= fit1$par[12],
                              stage3 = fit1$par[13], stage4 = fit1$par[14]),
               var = solve(-fit1$hessian)[11:14, 11:14],
               loglik = NA,
               score = NA,
               iter = NA,
               linear.predictors = NA,
               residuals = NA,
               means = NA,
               concordance = NA,
               method = "BFGS",
               n = length(delta),
               nevent = sum(delta),
               terms = NA,
               assign = NA,
               wald.test = coxph.wtest(solve(-fit1$hessian)[11:14, 11:14], fit1$par[11:14])$test,
               y = NA,
               formula = NA,
               call = cox$call
               )

class(out1) <- 'coxph'
xtable::xtable(out1, digits = 4)
out1
print((round(exp(out1$coefficients), 4) - round(exp(cox$coefficients), 4))/round(exp(cox$coefficients), 4) *100)

for (i in 1:4){
  print((out1$coefficients[i] + c(-1.96, + 1.96)*sqrt(out1$var[i, i])), digits = 4)
}
library(spsurv)
fit2 <- spbp(Surv(dat$time, dat$delta) ~ age + factor(stage),
             data = dat, approach = 'mle', hessian = TRUE)

out2 <- list(coefficients=c(age = fit2$par[1], stage2= fit2$par[2],
                              stage3 = fit2$par[3], stage4 = fit2$par[4]),
               var = solve(-fit2$hessian)[1:4, 1:4],
               loglik = NA,
               score = NA,
               iter = NA,
               linear.predictors = NA,
               residuals = NA,
               means = NA,
               concordance = NA,
               method = "BFGS",
               n = length(delta),
               nevent = sum(delta),
               terms = NA,
               assign = NA,
               wald.test = coxph.wtest(solve(-fit2$hessian)[1:4, 1:4], fit2$par[1:4])$test,
               y = NA,
               formula = NA,
               call = cox$call
)
class(out2) <- "coxph"
xtable::xtable(out2, digits = 4)
print(out2)
print((round(exp(out2$coefficients), 4) - round(exp(cox$coefficients), 4))/round(exp(cox$coefficients), 4) *100, digits = 5)

for (i in 1:4){
  print((out2$coefficients[i] + c(-1.96, + 1.96)*sqrt(out2$var[i, i])), digits = 4)
}

################################## BAYESIAN ESTIMATION ##################################
 rm(list = ls()); graphics.off(); library(survival);
 set.seed(1)
 options(scipen = 9999, OutDec = ".", digits = 4)
 library(KMsurv)
 data("larynx")
 dat <- larynx; rm(larynx)

dat$age <- as.numeric(scale(dat$age, scale = F))
dat$stage2 = as.numeric(scale(as.numeric(dat$stage == 2), scale = F))
dat$stage3 = as.numeric(scale(as.numeric(dat$stage == 3), scale = F))
dat$stage4 = as.numeric(scale(as.numeric(dat$stage == 4), scale = F))

cox <- coxph(Surv(dat$time, dat$delta) ~ age + stage2 + stage3 + stage4,
             data = dat)
source('fun.R')

# prior
a_gammak = .01
b_gammak = .01

a_gammak / b_gammak
## [1] 1
a_gammak / (b_gammak)^2
## [1] 100
m_beta <- rep(0, q)
S_beta <- diag(100, q, q)

## Gibbs sampling
t <- 1; it <- 2000
par_samp <- matrix(NA, it, m + q); colnames(par_samp) <- c(paste0("gamma", 1:m), paste0("beta", 1:q))
par <- c(rep(1, m), rep(0, q)); start_time <- Sys.time()
while(t < it + 1){
  for(k in 1:m){
    samp <- HI::arms(.5, myldens = function(y) l_gammak(y, k = k,
                                                        par = par,m = m, Z = Z, b = b, B = B, delta = delta),
                     indFunc = function(x){(x > 0) * (x < 1)}, n.sample = 1 )

    par_samp[t, k] <- -log(samp)
    par[k] <- par_samp[t, k]
  }
  for(j in (m + 1):(m + q)){
    samp <- HI::arms( .5, myldens = function(y) l_betaj(y, j = j,
                                                        par = par, m = m, Z = Z, b = b, B = B, delta = delta),
                      indFunc = function(x){(x > 0) * (x <1)}, n.sample = 1 )

    par_samp[t, j] <- log(samp) - log(1 - samp)
    par[j] <- par_samp[t, j]
  }
  cat( "\f", "Iteration", t, " out of ", it,". Time spent ",  Sys.time() - start_time)

  t <- t + 1
}
library(coda)
samp1 <- mcmc(par_samp[((it/2)+1):it,11:14]) # 16.46 sec
coda::traceplot(samp1[,1], main = '', bty = 'n')
coda::traceplot(samp1[,2], main = '', bty = 'n')
coda::traceplot(samp1[,3], main = '', bty = 'n')
coda::traceplot(samp1[,4], main = '', bty = 'n')

densplot(samp1[,1], main = '', bty = 'n')
densplot(samp1[,2], main = '', bty = 'n')
densplot(samp1[,3], main = '', bty = 'n')
densplot(samp1[,4], main = '', bty = 'n')

# exp_samp1 <- exp(samp1)
# traceplot(exp_samp1[,1], main = '', bty = 'n')
# traceplot(exp_samp1[,2], main = '', bty = 'n')
# traceplot(exp_samp1[,3], main = '', bty = 'n')
# traceplot(exp_samp1[,4], main = '', bty = 'n')
#
# densplot(exp_samp1[,1], main = '', bty = 'n')
# densplot(exp_samp1[,2], main = '', bty = 'n')
# densplot(exp_samp1[,3], main = '', bty = 'n')
# densplot(exp_samp1[,4], main = '', bty = 'n')

library(coda)
xtable::xtable(cbind(Median = apply(samp1, 2, median), summary(samp1)[[1]],
                     HPDL = coda::HPDinterval(samp1)[,1], HPDU = coda::HPDinterval(samp1)[,2])[,c(1,2,3,6,7)], digits = 4)
xtable::xtable(cbind(Median = apply(exp_samp1, 2, median), summary(exp_samp1)[[1]],
                     HPDL = coda::HPDinterval(exp_samp1)[,1], HPDU = coda::HPDinterval(exp_samp1)[,2])[,c(1,2,3,6,7)], digits = 4)

# print((round(apply(samp1, 2, mean), 4) - round(cox$coefficients, 4))/round(cox$coefficients, 4) *100, digits = 5)

round(exp(round(apply(samp1, 2, mean), 4)),4)

pdf('tracebeta1_arms_larynx.pdf')
coda::traceplot(samp1[,1], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
dev.off()
pdf('tracebeta2_arms_larynx.pdf')
coda::traceplot(samp1[,2], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
dev.off()
pdf('tracebeta3_arms_larynx.pdf')
coda::traceplot(samp1[,3], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
dev.off()
pdf('tracebeta4_arms_larynx.pdf')
coda::traceplot(samp1[,4], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
dev.off()

pdf('densbeta1_arms_larynx.pdf')
coda::densplot(samp1[,1], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
dev.off()
pdf('densbeta2_arms_larynx.pdf')
coda::densplot(samp1[,2], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
dev.off()
pdf('densbeta3_arms_larynx.pdf')
coda::densplot(samp1[,3], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
dev.off()
pdf('densbeta4_arms_larynx.pdf')
coda::densplot(samp1[,4], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
dev.off()

# pdf('traceHR1_arms_larynx.pdf')
# coda::traceplot(exp_samp1[,1], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
# dev.off()
# pdf('traceHR2_arms_larynx.pdf')
# coda::traceplot(exp_samp1[,2], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
# dev.off()
# pdf('traceHR3_arms_larynx.pdf')
# coda::traceplot(exp_samp1[,3], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
# dev.off()
# pdf('traceHR4_arms_larynx.pdf')
# coda::traceplot(exp_samp1[,4], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
# dev.off()

# pdf('densHR1_arms_larynx.pdf')
# coda::densplot(exp_samp1[,1], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
# dev.off()
# pdf('densHR2_arms_larynx.pdf')
# coda::densplot(exp_samp1[,2], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
# dev.off()
# pdf('densHR3_arms_larynx.pdf')
# coda::densplot(exp_samp1[,3], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
# dev.off()
# pdf('densHR4_arms_larynx.pdf')
# coda::densplot(exp_samp1[,4], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
# dev.off()

require(spsurv)
fitbe <- spbp(Surv(dat$time, dat$delta) ~ age + factor(stage),
     data = dat, approach = 'bayes', chain = 1)

library(coda)
samp2 <- mcmc(rstan::extract(fitbe, pars = "beta")$beta) # 13.93 sec.
xtable::xtable(cbind(Median = apply(samp2, 2, median), summary(samp2)[[1]],
                     HPDL = coda::HPDinterval(samp2)[,1], HPDU = coda::HPDinterval(samp2)[,2])[,c(1,2,3,6,7)], digits = 4)

print((round(apply(samp2, 2, mean), 4) - round(cox$coefficients, 4))/round(cox$coefficients, 4) *100, digits = 4)

round(exp(round(apply(samp2, 2, mean), 4)),4)

pdf('tracebeta1_spbp_larynx.pdf')
coda::traceplot(samp2[,1], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
dev.off()
pdf('tracebeta2_spbp_larynx.pdf')
coda::traceplot(samp2[,2], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
dev.off()
pdf('tracebeta3_spbp_larynx.pdf')
coda::traceplot(samp2[,3], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
dev.off()
pdf('tracebeta4_spbp_larynx.pdf')
coda::traceplot(samp2[,4], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
dev.off()

pdf('densbeta1_spbp_larynx.pdf')
coda::densplot(samp2[,1], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
dev.off()
pdf('densbeta2_spbp_larynx.pdf')
coda::densplot(samp2[,2], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
dev.off()
pdf('densbeta3_spbp_larynx.pdf')
coda::densplot(samp2[,3], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
dev.off()
pdf('densbeta4_spbp_larynx.pdf')
coda::densplot(samp2[,4], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
dev.off()
print(fitbe)
