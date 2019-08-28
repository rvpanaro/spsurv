// include Functions block.
#include /chunks/loglikbp.stan

// Data block (important).
data{
  // setting observed data:
  int<lower=1> n;
  int<lower=1> m;
  int<lower=1> q;
  real<lower= m> tau;
  int<lower=0, upper=1> approach;
  int<lower=0, upper=1> null;
  int<lower=0, upper=2> M;
  vector<lower=0, upper = 1>[n] status;
  vector<lower=0>[n] time;

  // observed arrays:
  matrix[n,q] X;
  matrix[n,m] b;
  matrix[n,m] B;

  // setting hyperparameters:
  real<lower=0> shape_gamma;
  real<lower=0> rate_gamma;
  real mean_beta;
  real sd_beta;
}

// Parametes block (important).
parameters{
vector[q] beta;
vector<lower=0>[m] gamma;
}
// Model block (important).
model{
  vector[n] loglik;
    if(M == 0){
      loglik = loglikpo(beta, gamma, status, X, b, B, null);
    }
    else if( M == 1){
      loglik = loglikph(beta, gamma, status, X, b, B, null);
    }
    else{
      loglik = loglikaft(time, beta, gamma, status, X, b, B, tau, null);
    }

    target += sum(loglik);

  if( approach == 1){
	  beta ~ normal(mean_beta,sd_beta);
	  gamma ~ gamma(shape_gamma, rate_gamma);
  }
}
// Final line empty to avoid warnings.
