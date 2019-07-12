## Proportional Odds Example according to Bennet (1983)
rm(list = ls(all = T)); library(survival); data(veteran)

## Only 97 patients with no prior therapy are considered
veteran <- veteran[veteran$prior == 0,]

# Patients with high (over 50) and low performance status (PS)
veteran$PS <- cut(veteran$karno,c(0,50, 100), labels = c("Low", "High"))

# Cell type
veteran$celltype <- factor(veteran$celltype, levels = c("large", "adeno", "smallcell", "squamous"))

# Kaplan Meier
kme_PS <- survfit(formula = Surv(time, status) ~ PS, data = veteran);

# pdf('kmeodds.pdf'); dev.off()
plot(kme_PS, bty = 'n', xlab = "time",lwd = 3, cex.lab = 1.7, cex.axis = 1.7)
kme_celltype <- survfit(formula = Surv(time, status) ~ celltype, data = veteran);
plot(kme_celltype, xlab = "time",lwd = 3, cex.lab = 1.7, cex.axis = 1.7, bty = "n")

plot((1-kme_celltype[1]$surv)/kme_celltype[1]$surv, pch = 19, bty = 'n')
points((1-kme_celltype[2]$surv)/kme_celltype[2]$surv, pch = 2)

plot((1-kme_celltype[1]$surv)/kme_celltype[1]$surv, pch = 19, bty = 'n')
points(log((1-kme_celltype[3]$surv)/kme_celltype[3]$surv), pch = 3)

plot((1-kme_celltype[1]$surv)/kme_celltype[1]$surv, pch = 19, bty = 'n')
points((1-kme_celltype[4]$surv)/kme_celltype[4]$surv, pch = 4)
###

### Making use of spsurv package
library(spsurv)
# devtools::install_github('rvpanaro/spsurv')
options(scipen=9999)

library(timereg)

# Transformação em dados intervalares
veteran$time_start = veteran$time
veteran$time_end = if_else(veteran$status == 1, veteran$time, Inf)

fit <- prop.odds(Event(time, status) ~ karno + factor(celltype), data = veteran)
summary(fit)

plot(predict.timereg(fit, Z = c(60,0,0,0))$time, predict.timereg(fit, Z = c(0,0,0,0))$S0, bty = 'n', lwd = 4, cex.lab = 1.7, cex.axis = 1.7, xlab = 'time', ylab = '',type = 's')
lines(predict.timereg(fit, Z = c(60,1,0,0))$time, predict.timereg(fit, Z = c(0,1,0,0))$S0, lwd = 4, type = 's')
lines(predict.timereg(fit, Z = c(60,0,1,0))$time, predict.timereg(fit, Z = c(0,0,1,0))$S0, lwd = 4, type = 's')
lines(predict.timereg(fit, Z = c(60,0,0,1))$time, predict.timereg(fit, Z = c(0,0,0,1))$S0, lwd = 4, type = 's')

for (i in 1:4){
  print(exp(fit$gamma[[i]] + c(-1.96, + 1.96)*sqrt(fit$var.gamma[i, i])), digits = 4)
}

#1) Parameter estimates (and standard errors) for the full model
# were - 0.053 (0.010) for PS;
#1.31 (0.55) for adeno vs. large cell type;
#1.38 (0.52) for small vs. large cell type
# and - 0.18 (0.59)
# for squamous vs. large cell type.

# dat$celltype <- relevel(dat$celltype, "large")
fitspbp<- spsurv::spbp(Surv(time, status) ~ karno + celltype, dat = veteran,
                        approach = "mle", model = 'po', hessian = T)
fitspbp
S <- function(y, beta, gamma, tau, z){
  m <- length(gamma)
  k <- 1:m
  B <- sapply(k, function(k) pbeta(y/tau, k, m - k + 1))
  BP <- apply(B, 1, function(x){gamma %*% x})
  theta <- exp(t(z) %*% beta)
  S <- 1/(1 + as.vector(theta) * BP)
  return(S)
}
beta = fit1$par[1:4]
gamma = fit1$par[5:14]
curve(S(x, beta = beta, gamma = gamma, tau = max(veteran$time), z = c(0,0,0,0)), lwd = 4, col = "darkblue", add = T)
curve(S(x, beta = beta, gamma = gamma, tau = max(veteran$time), z = c(0,1,0,0)), lwd = 4, col = "orange", add = T)
curve(S(x, beta = beta, gamma = gamma, tau = max(veteran$time), z = c(0,0,1,0)), lwd = 4, col = "darkblue", add = T, lty = 2)
curve(S(x, beta = beta, gamma = gamma, tau = max(veteran$time), z = c(0,0,0,1)), lwd = 4, col = "darkblue", add = T, lty = 3)

library(survival)
out1 <- list(coefficients=c(karno = fit1$par[1], celltypeadeno= fit1$par[2],
                            celltypesmallcell = fit1$par[3], celltypesquamous = fit1$par[4]),
             var = solve(-fit1$hessian)[1:4, 1:4],
             loglik = NA,
             score = NA,
             iter = NA,
             linear.predictors = NA,
             residuals = NA,
             means = NA,
             concordance = NA,
             method = "BFGS",
             n = length(veteran$status),
             nevent = sum(veteran$status),
             terms = NA,
             assign = NA,
             wald.test = coxph.wtest(solve(-fit1$hessian)[1:4, 1:4], fit1$par[1:4])$test,
             y = NA,
             formula = NA,
             call = "spsurv::spbp(Surv(time, status) ~ karno + celltype, dat = dat,
                        approach = 'mle', model = 'po', hessian = T)"
)

class(out1) <- 'coxph'
xtable::xtable(out1, digits = 4)
print((round(exp(out1$coefficients), 4) - round(exp(fit$gamma), 4))/round(exp(fit$gamma), 4) *100)
out1
for (i in 1:4){
  print(exp(out1$coefficients[i] + c(-1.96, + 1.96)*sqrt(out1$var[i, i])), digits = 4)
}

fit1be <- spsurv::spbp(Surv(time, status) ~ karno + celltype, dat = veteran,
                       approach = "bayes", model = 'po', chain  = 1)

print(fit1be, par = 'beta')
# rstan::traceplot(fit1be, par = 'beta')

library(coda)
samp2 <- mcmc(rstan::extract(fit1be, pars = "beta")$beta) # 13.93 sec.
xtable::xtable(cbind(Median = apply(samp2, 2, median), summary(samp2)[[1]],
                     HPDL = coda::HPDinterval(samp2)[,1], HPDU = coda::HPDinterval(samp2)[,2])[,c(1,2,3,6,7)], digits = 4)

# pdf('tracebeta1_spbp_lung.pdf')
coda::traceplot(samp2[,1], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
# dev.off()
# pdf('tracebeta2_spbp_lung.pdf')
coda::traceplot(samp2[,2], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
# dev.off()
# pdf('tracebeta3_spbp_lung.pdf')
coda::traceplot(samp2[,3], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
# dev.off()
# pdf('tracebeta4_spbp_lung.pdf')
coda::traceplot(samp2[,4], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
# dev.off()

# pdf('densbeta1_spbp_lung.pdf')
coda::densplot(samp2[,1], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
# dev.off()
# pdf('densbeta2_spbp_lung.pdf')
coda::densplot(samp2[,2], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
# dev.off()
# pdf('densbeta3_spbp_lung.pdf')
coda::densplot(samp2[,3], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
# dev.off()
# pdf('densbeta4_spbp_lung.pdf')
coda::densplot(samp2[,4], bty = 'n', cex.lab = 2.5, cex.axis = 2.5)
# dev.off()

est <- c(-0.053, 1.31, 1.38, -0.18)
cbind(est, fit1mle$par[1:4])
cat( "Percentage Error%\n",
     ((fit1mle$par[1:4]-est)/abs(est)) * 100)

xtable(cbind(round(cbind(est,  fit1mle$par[1:4]),4), ((fit1mle$par[1:4]-est)/abs(est)) * 100), digits = 4)

print(fit1be, digits_summary = 4)

f1be <- c(-0.0536, 1.3057, 1.3634, -0.1617)
xtable(cbind(round(cbind(est,  f1be),4), ((f1be-est)/abs(est)) * 100), digits = 4)

cat( "Percentage Error%\n",
     ((f1be-est)/abs(est)) * 100)

###########################################################
##### ARMS within gibbs sampling ##########################
###########################################################
rm(list = ls())
graphics.off();
set.seed(1)
options(scipen = 9999, OutDec = ".", digits = 4)
library("survival")
data("veteran")
veteran <- veteran[veteran$prior == 0,]
veteran$celltype <- factor(veteran$celltype, levels = c("large", "adeno", "smallcell", "squamous"))
cox <- coxph(Surv(time, status) ~ karno + celltype,
             data = veteran)
library(timereg)

fit <- prop.odds(Event(time, status) ~ karno + celltype, data = veteran)
fitspbp <- spsurv::spbp(Surv(time, status) ~ karno + celltype, dat = veteran,
                       approach = "mle", model = 'po', hessian = T)

setwd("~/Documents/spsurv/inst/PO")
m =10
source('fun.R')
fit1 <- optim(c(rep(1, m), rep(0, q)), fn = loglik, method = "BFGS", control = list(fnscale = -1, maxit = 2000), m = m,
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
print(out1)

print((round(exp(out1$coefficients), 4) - round(exp(fit$gamma), 4))/round(exp(fit$gamma), 4) *100)

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

xtable::xtable(cbind(Median = apply(samp1, 2, median), summary(samp1)[[1]],
                     HPDL = coda::HPDinterval(samp1)[,1], HPDU = coda::HPDinterval(samp1)[,2])[,c(1,2,3,6,7)], digits = 4)
coef <- exp(c(-0.0533, 1.3307, 1.4123, -0.1556))
print((coef-round(exp(apply(samp1, 2, mean)), 4))/round(exp(apply(samp1, 2, mean)), 4)  *100, digits = 4)
