// include Functions block.
#include /include/loglikbp.stan

// Data block (important).
data{
  // setting observed data:
  int<lower=1> n;
  int<lower=1> m;
  int<lower=1> q;
  real tau;
  int<lower=0, upper=1> approach;
  int<lower=0, upper=3> dist;
  int<lower=0, upper=1> null;
  int<lower=0, upper=2> M;
  vector<lower=0, upper = 1>[n] status;
  int<lower=0> id [n];
  vector<lower=0>[n] z; // for coding purposes
  vector<lower=0>[n] time;

  // observed arrays:
  matrix[n,q] X;
  matrix[n,m] b;
  matrix[n,m] B;

  // setting hyperparameters:
  real<lower=0> hyper_gamma;
  real mean_beta;
  real<lower=0> sd_beta;
  real mean_nu;
  real<lower=0> sd_nu;
}

// Parametes block (important).
parameters{
  vector[q] beta;
  vector[m] nu;
}

transformed parameters{
    vector[n] log_lik;

    if(null == 1){
      for(i in 1:n){
        log_lik = loglik_null(beta, nu, status, X, b, B, M, dist, id, z);
      }
    }
    else{
      if(M == 0){
          log_lik = loglik_po(beta, nu, status, X, b, B, dist, id, z);
      }
      else if( M == 1){
          log_lik = loglik_ph(beta, nu, status, X, b, B, dist, id, z);
      }
      else{
          log_lik = loglik_aft(time, beta, nu, status, X, b, B, dist, id, z);
      }
    }
}

// Model block (important).
model{

  if( approach == 1){
	  beta ~ normal(mean_beta, sd_beta);
	  nu ~ normal(mean_nu, sd_nu);
   }

   target += sum(log_lik);
}


// Final line empty to avoid warnings.

