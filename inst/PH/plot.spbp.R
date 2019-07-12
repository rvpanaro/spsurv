require(survival)
H <- function(y, beta, gamma, tau, z){
  m <- length(gamma)
  k <- 1:m
  B <- sapply(k, function(k) pbeta(y/tau, k, m - k + 1))
  H0 <- apply(B, 1, function(x){gamma %*% x})
  theta <- exp(t(z) %*% beta)
  H <- as.vector(theta) * H0
  return(H)
}

Ht<-basehaz(cox, centered = T)
times <-Ht$time
H0 <-Ht$hazard ## Breslow non-parametric Estimate
S0 <- exp(-H0)

beta <- fit2$par[c(1,2,3,4)]
gamma <- fit2$par[-c(1,2,3,4)]

plot(Ht$time, S0, type = 'l', bty = 'n')

curve(exp(-H(x, beta = beta, gamma = gamma, tau = tau, z = c(0,0,0,0))), lwd = 2, lty = 2, add = T)

new <- data.frame(age = c(77-attr(dat$age, 'scaled:center'),
                          77-attr(dat$age, 'scaled:center'),
                          77-attr(dat$age, 'scaled:center'),
                          77-attr(dat$age, 'scaled:center')),
                  stage2 = c(0-attr(dat$stage2, 'scaled:center'),
                             1-attr(dat$stage2, 'scaled:center'),
                             0-attr(dat$stage2, 'scaled:center'),
                             0-attr(dat$stage2, 'scaled:center')),
                  stage3 = c(0-attr(dat$stage3, 'scaled:center'),
                             0-attr(dat$stage3, 'scaled:center'),
                             1-attr(dat$stage3, 'scaled:center'),
                             0-attr(dat$stage3, 'scaled:center')),
                  stage4 = c(0-attr(dat$stage4, 'scaled:center'),
                             0-attr(dat$stage4, 'scaled:center'),
                             0-attr(dat$stage4, 'scaled:center'),
                             1-attr(dat$stage4, 'scaled:center')))

st_cox <- survfit(cox, newdata = new)

pdf('surv_mle_larynx.pdf')
plot(st_cox, bty = 'n', lwd = 4, xlab = 'time', cex.axis = 1.7, cex.lab = 1.7, col = c(1,2,3,4))
curve(exp(-H(x, beta = beta, gamma = gamma, tau = tau, z = as.numeric(new[1,]))), lwd = 4, lty = 4, add = T, col = 1)
curve(exp(-H(x, beta = beta, gamma = gamma, tau = tau, z = as.numeric(new[2,]))), add = T, lwd = 4, lty = 4, col = 2)
curve(exp(-H(x, beta = beta, gamma = gamma, tau = tau, z = as.numeric(new[3,]))), add = T, lwd = 4, lty = 4, col = 3)
curve(exp(-H(x, beta = beta, gamma = gamma, tau = tau, z = as.numeric(new[4,]))), add = T, lwd = 4, lty = 4, col = 4)
dev.off()

### Bayesian
set.seed(1)
require(spsurv)
fitbe <- spbp(Surv(dat$time, dat$delta) ~ age + factor(stage),
              data = dat, approach = 'bayes', chain = 1)
library(coda)
library(rstan)
par_samp <- cbind(rstan::extract(fitbe, pars = "beta")$beta,
                  rstan::extract(fitbe, pars = "gamma")$gamma)# 13.93 sec.

H <- function(par, Z, B){
  q <- length(par) - m
  # print(q); print(m)
  beta <- matrix(par[1:q], ncol = 1)
  gamma <- matrix(par[(q + 1):(m + q)], ncol = 1)
  gamma <- abs(gamma)
  # print(dim(Z)); print(dim(beta))
  eta <- Z %*% beta
  # h <-  b %*% gamma * exp (eta)
  H <- B %*% gamma * exp (eta)
}

H_bern_e0 <- apply(apply(par_samp, 1, H, Z = matrix( c(rep(12.388889,90), rep(-0.1888889,90), rep(-0.3,90), rep(-0.1444444,90) )
                                                     , nrow = 90, ncol = 4), B =B), 1, mean)
H_bern_e1 <- apply(apply(par_samp, 1, H, Z = matrix( c(rep(12.388889,90), rep(0.8111111,90), rep(-0.3,90), rep(-0.1444444,90) )
                                                     , nrow = 90, ncol = 4), B =B), 1, mean)
H_bern_e2 <- apply(apply(par_samp, 1, H, Z = matrix( c(rep(12.388889,90), rep(-0.1888889,90), rep(0.7,90), rep(-0.1444444,90) )
                                                     , nrow = 90, ncol = 4), B =B), 1, mean)
H_bern_e3 <- apply(apply(par_samp, 1, H, Z = matrix( c(rep(12.388889,90), rep(-0.1888889,90), rep(-0.3,90), rep(0.8555556,90) )
                                                      , nrow = 90, ncol = 4), B =B), 1, mean)
pdf('surv_bayes_larynx.pdf')
plot(st_cox, bty = 'n', lwd = 4, xlab = 'time', cex.axis = 1.7, cex.lab = 1.7, col = c(1,2,3,4))
lines(sort(larynx$time), exp(-H_bern_e0[order(larynx$time)]), col = 1, lwd = 4, lty = 4)
lines(sort(larynx$time), exp(-H_bern_e1[order(larynx$time)]), col = 2, lwd = 4, lty = 4)
lines(sort(larynx$time), exp(-H_bern_e2[order(larynx$time)]), col = 3, lwd = 4, lty = 4)
lines(sort(larynx$time), exp(-H_bern_e3[order(larynx$time)]), col = 4, lwd = 4, lty = 4)
dev.off()
