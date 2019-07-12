set.seed(1)
source("ph_gendata.R")

# Generated sample size
n <- 200
#Coefficients
beta_init <- c(-1,1,3,-2) # beginning with beta0
#Covariates
Z <- cbind(rbinom(n = n, size = 1, prob = .5),
           rnorm(n = n, 0, 1),
           runif(n = n, -2, 2))

dat <- exp_surv_gen(n, beta_init, Z)
saveRDS(data, 'simdata.R')

require(survival); require(xtable); require(spsurv)
cox <- coxph(Surv(time, status) ~ X1 + X2 + X3, data = dat)
xtable(cox, digits = 4)
print((round(cox$coefficients, 4) - round(beta_init[-1], 4))/round(beta_init[-1], 4) *100, digits = 4)

source('fun.R')

# Classical
fit1 <- optim(c(rep(1, m), rep(0, q)), fn = loglik, method = "BFGS", control = list(fnscale=-1), m = m,
              delta = dat$status, Z = Z, b = b, B = B, hessian = TRUE)

out1 <- list(coefficients=c(x1 = fit1$par[16], x2= fit1$par[17],
                            x3 = fit1$par[18]),
             var = solve(-fit1$hessian)[16:18, 16:18],
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
             wald.test = coxph.wtest(solve(-fit1$hessian)[16:18, 16:18], fit1$par[16:18])$test,
             y = NA,
             formula = NA,
             call = cox$call
)

class(out1) <- 'coxph'
xtable::xtable(out1, digits = 4)
print((round(out1$coefficients, 4) - round(cox$coefficients, 4))/round(cox$coefficients, 4) *100)

fitmle <- spbp(Surv(time, status) ~ X1 + X2 +X3, data = dat,
            approach = 'mle', hessian = T)

out2 <- list(coefficients=fitmle$par[1:3],
               var = solve(-fitmle$hessian)[1:3, 1:3],
               loglik = NA,
               score = NA,
               iter = NA,
               linear.predictors = NA,
               residuals = NA,
               means = NA,
               concordance = NA,
               method = "BFGS",
               n = length(dat$status),
               nevent = sum(dat$status),
               terms = NA,
               assign = NA,
               wald.test = coxph.wtest(solve(-fitmle$hessian)[1:3, 1:3], fitmle$par[1:3])$test,
               y = NA,
               formula = NA,
               call = cox$call
)

class(out2) <- 'coxph'
xtable::xtable(out2, digits = 4)
print((round(fitmle$par[1:3], 4) - round(beta_init[-1], 4))/round(beta_init[-1], 4) *100, digits = 5)
##################### Bayesian
str(data)

Z <- matrix(model.matrix(cox), nrow = 200); colnames(Z) <- names(cox$coefficients)
n <- nrow(Z)
q <- ncol(Z)
m <- ceiling(sqrt(n)); cat("m = ", m)
tau <- max(data$time); cat("tau = ", tau)
delta = data$delta

# Polynomial basis
k <- 1:m
b <- matrix(NA, n, m )
B <- matrix(NA, n, m )
y <- data$time/tau

for (i in 1:n){
  for(k in 1:m){
    b[i,k] <- dbeta(y[i], k, m - k + 1) / tau
    B[i,k] <- pbeta(y[i], k, m - k + 1)
  }
}

loglik <- function(par, ## variables
                   m, delta, Z, b, B){ ## data

  q <- length(par) - m
  gamma <- matrix(par[1:m], ncol = 1)
  gamma <- abs(gamma)
  beta <- matrix(par[(m + 1):(m + q)], ncol = 1)
  eta <- Z %*% beta

  h <-  b %*% gamma
  H <-  B %*% gamma
  res <- sum(delta * log(h * exp(eta)) - H * exp(eta))

  return(res)
}

# prior
a_gammak = .01
b_gammak = .01

a_gammak / b_gammak
## [1] 1
a_gammak / (b_gammak)^2
## [1] 100
m_beta <- rep(0, q)
S_beta <- diag(100, q, q)

# Bayesian
## Conditionals
l_gammak <- function(pi_k, k, ## variables
                     par, m, delta, Z, b, B){  ## data

  gammak <- -log(pi_k); res <- NULL

  f <- function(par) {
    loglik( par = par, m = m, Z = Z, delta = delta, b = b, B = B )  +
      + (a_gammak - 1) * log(par[k]) - b_gammak * par[k]
  }

  # gammak = 2; par = c(1,1,1,1,1,1,1,1,11,1,0,0,10,10); k = 3
  # gammak <- c(1,2,3)
  for ( i in 1:length(gammak)){
    par_aux <- par
    par_aux[k] <- gammak[i]
    res[i] <- f(par_aux)
  }

  return( res - log(pi_k))
}

# l_gammak( c(.3,.1,.14,.81,.91), par = c(1,1,1,1,1,1,1,1,11,1,0,0,10,100), k = 8, m = m, Z = Z, b = b, B = B, delta = delta)
l_betaj <- function(pi_j, j, ## variables
                    par, m, delta, Z, b, B){ ## data

  # pi_j <- c(.30,.1,.14,.81,.91)
  # j <- 12
  q <- j - m

  betaj <- log(pi_j) - log(1 - pi_j); res <- NULL

  f <- function(par) {

    loglik( par = par, m = m, Z = Z, delta = delta, b = b, B = B ) -
      ( 1 / 2 ) * (par[j] - m_beta[q]) ^ 2   * 1/S_beta[q, q]
  }
  # gammak <- c(1,2,3)
  for ( i in 1:length(betaj)){
    par_aux <- par
    par_aux[j] <- betaj[i]
    res[i] <- f(par_aux)
  }
  return(res - log(pi_j) - log(1 - pi_j))
}

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
fitmle2 <- optim(c(rep(1, m), rep(0, q)), fn = loglik, method = "BFGS", control = list(fnscale=-1), m = m,
             delta = delta, Z = Z, b = b, B = B, hessian = TRUE)

require(coda)
samp2 <- mcmc(par_samp[,16:18])

xtable::xtable(cbind(Median = apply(samp2, 2, median), summary(samp2)[[1]],
                   HPDL = coda::HPDinterval(samp2)[,1], HPDU = coda::HPDinterval(samp2)[,2])[,c(1,2,3,6,7)], digits = 4)

fitbe <- spbp(Surv(time, status) ~ X1 + X2 +X3, data = data, chain = 1,
              priors = list(shape_gamma = 0.1, rate_gamma = 0.1, mean_beta = 0,
                            sd_beta = 10))
library(rstan)
samp1 <- mcmc(extract(fitbe)$beta)
xtable::xtable(cbind(Median = apply(samp1, 2, median), summary(samp1)[[1]],
                     HPDL = coda::HPDinterval(samp1)[,1], HPDU = coda::HPDinterval(samp1)[,2])[,c(1,2,3,6,7)], digits = 4)
