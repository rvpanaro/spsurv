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
  int<lower=0> priordist [2];
  real priorpars [4];

  // beta (vector)
  int<lower=0> priordist_beta [q];
  real location_beta [q];
  real<lower=0> scale_beta [q];

  // Standard quantities
  vector<lower=0>[q] std; // feature standard deviatons
  vector[q] means; // feature means
}

// Parametes block (important).
parameters{
  vector[q] beta; // feature effect
  vector<lower=0>[m] gamma; //BP basis effect
}

transformed parameters{
    vector[n] log_lik; // log likelihood
    vector[m] nu; // exp of the BP basis effect

    vector[q] beta_std; // standard beta
    vector<lower=0>[m] gamma_std; // standard BP gamma

    beta_std = beta ./ std;

    if(M == 2){
        gamma_std = gamma * exp(sum(beta .* means ./ std));
    }
    else{
        gamma_std = gamma * exp(-sum(beta .* means ./ std));
    }

    nu = log(gamma);

    if(null == 1){
      for(i in 1:n){
        log_lik = loglik_null(beta, gamma, status, X, b, B, M, dist, id, z);
      }
    }
    else{
      if(M == 0){
          log_lik = loglik_po(beta, gamma, status, X, b, B, dist, id, z);
      }
      else if( M == 1){
          log_lik = loglik_ph(beta, gamma, status, X, b, B, dist, id, z);
      }
      else{
          log_lik = loglik_aft(time, beta, gamma, status, X, b, B, dist, id, z);
      }
    }
}

// Model block (important).
model{

if(approach == 1){ // priors
////// Beta Prior
    for(i in 1:q){
      // print("location = ", location_beta[i], " scale = ", scale_beta[i]);

      if(priordist_beta[i] == 0){
        beta ~ normal(location_beta[i], scale_beta[i]);
      }
      else{
        beta ~ cauchy(location_beta[i], scale_beta[i]);
      }
    }
////// Gamma Prior
    // print("h1 = ", priorpars[1], " h2 = ", priorpars[2]);
    if(priordist[1] == 1){
      gamma ~ gamma(priorpars[1], priorpars[2]);
    }
    else if(priordist[1] == 2){
      gamma ~ inv_gamma(priorpars[1], priorpars[2]);
    }
    else{
      gamma ~ lognormal(priorpars[1], priorpars[2]);
    }
}

   target += sum(log_lik);
}

// Final line empty to avoid warnings.

