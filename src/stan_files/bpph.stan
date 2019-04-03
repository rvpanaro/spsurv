// include Functions block.
#include /chunks/loglikbp.stan

// Data block (important).
data{
  // setting observed data:
  int<lower=1> n;
  int<lower=1> m;
  int<lower=1> q;
  int<lower=0, upper=1> approach;
  vector<lower=0, upper = 1>[n] status;

  // observed arrays:
  matrix[n,q] Z;
  matrix[n,m] b;
  matrix[n,m] B;

  // setting hyperparameters:
  real<lower=0> a_gamma;
  real<lower=0> b_gamma;
  real m_beta;
  real S_beta;
}

// Parametes block (important).
parameters{
vector[q] beta;
vector<lower=0>[m] gamma;
}
// Model block (important).
model{
  vector[n] loglik;
  loglik = loglikbp(beta, gamma, status, Z, b, B);
  target += sum(loglik);

  if( approach == 1){
	beta ~ normal(m_beta,S_beta);
	gamma ~ gamma(a_gamma, b_gamma);
  }
}
// Final line empty to avoid warnings.
