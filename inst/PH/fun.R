# Data
Z <- matrix(model.matrix(cox), nrow = dim(dat)[1]); colnames(Z) <- names(cox$coefficients)
Z <-  scale(Z, scale = F)
n <- nrow(Z)
q <- ncol(Z)
m <- ceiling(sqrt(n)); cat("m = ", m)
## m =  10
tau <- max(dat$time); cat("tau = ", tau)
## tau =  10.7
delta = dat$delta
if(is.null(delta)){delta = dat$status}

# Polynomial basis
k <- 1:m
b <- matrix(NA, n, m )
B <- matrix(NA, n, m )
y <- dat$time/tau

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
