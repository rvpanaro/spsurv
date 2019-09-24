// include Functions block.
#include /include/loglikbp.stan

// Data block (important).
data{
  // setting observed data:
  int<lower=1> n;
  int<lower=1> m;
  int<lower=1> q;
  //int<lower=0> id[n]; // used for frailty models
  real tau;
  int<lower=0, upper=1> approach;
  int<lower=0, upper=3> dist;
  int<lower=0, upper=1> null;
  int<lower=0, upper=2> M;
  vector<lower=0, upper = 1>[n] status;
  int<lower=0> id [n];
  vector<lower=0>[n] time;

  // observed arrays:
  matrix[n,q] X;
  matrix[n,m] b;
  matrix[n,m] B;

  // setting hyperparameters:
  real<lower=0> shape_gamma;
  real<lower=0> rate_gamma;
  real mean_beta;
  real<lower=0> sd_beta;
}

// Parametes block (important).
parameters{
  vector[q] beta;
  vector[m] nu;
  vector[n] z;
}

transformed parameters{
   vector<lower=0>[n] exp_z = exp(z);
}

// Model block (important).
model{
  vector[n] loglik;
  if(null == 1){
      loglik = loglik_null(beta, nu, status, X, b, B, M, dist, id, z);
  }
  else{
    if(M == 0){
      loglik = loglik_po(beta, nu, status, X, b, B, dist, id, z);
    }
    else if( M == 1){
      loglik = loglik_ph(beta, nu, status, X, b, B, dist, id, z);
    }
    else{
      loglik = loglik_aft(time, beta, nu, status, X, b, B, dist, id, z);
    }
  }
    target += sum(loglik);

  if( approach == 1){
	  beta ~ normal(mean_beta, sd_beta);
	  nu ~ normal(mean_beta, sd_beta);

	  //frailty
	  if(dist == 1){
	    exp_z ~ gamma(0.01, 0.01); // gamma frailty
	  }
	  else if(dist == 2){
	    z ~ normal(0, 0.01); // gaussian frailty
	  }
	  else if(dist == 3){
	    z ~ student_t(30, 0, 0.01); // t-student frailty
	  }
	  else{
	    z ~ normal(0, 0.000001);
	  }
  }
}
// Final line empty to avoid warnings.
