// include Functions block.
#include /include/loglikbp.stan

// Data block (important).
data{
  // setting observed data:
  int<lower=1> n;
  int<lower=1> m;
  int<lower=1> q;
  real tau;
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
  vector[n] z;
  real<lower=0> kappa;
}

transformed parameters {
  real<lower=0> sigma;
  sigma = inv_sqrt(kappa);
}

// Model block (important).
model{

	  //priors
	  beta ~ normal(mean_beta, sd_beta);
	  nu ~ normal(mean_nu, sd_nu);
    kappa ~ gamma(hyper_gamma, hyper_gamma);

	  //frailty
	  if(dist == 1){
	    z ~ gamma(kappa, kappa); // gamma frailty
	  }
	  else if(dist == 2){
	    z ~ normal(0, sigma); // gaussian frailty
	  }
	  else if(dist == 3){
	    z ~ student_t(20, 0, sigma); // t-student frailty
	  }
	  else{
	    z ~ normal(0, 0.000001);
	  }
}

// Final line empty to avoid warnings.


