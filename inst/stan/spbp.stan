// include Functions block.
#include /include/loglikbp.stan

// Data block (important).
data{
  // observed data:
  int<lower=1> n;       //sample size
  int<lower=1> m;       // BP degree
  int<lower=1> p;       // number of features
  // real tau;       // right limit
  int<lower=0, upper=1> approach;       // 0 - freq or 1 - bayes
  int<lower=0, upper=2> rand;       // 0 - no effect, 1 - multiplicative or 2 - additive random effect
  int<lower=0, upper=1> null;       // 0 - null 1 - not null (set internally)
  int<lower=0, upper=2> M;       // 0 - po, 1 - ph or 2 - aft
  vector<lower=0, upper = 1>[n] status;       // censoring indicator ( 0 or 1 for each obs)
  int<lower=0> id [n];       // random effect identificator
  vector<lower=0>[n] z;       // random effect latent variable
  vector<lower=0>[n] time;       // response variable

  // observed arrays:
  matrix[n,p] X;       // features
  matrix[n,m] g;       // BP basis
  matrix[n,m] G;       // cumulative BP basis

  // prior distribution & hyperparameters:

  //// beta (vector)
  int<lower=0> priordist_beta [p];
  vector [p] location_beta;
  vector<lower=0> [p] scale_beta;

  //// gamma (vector)
  int<lower=0> priordist_gamma [m];
  vector<lower=0> [m] par1_gamma;
  vector<lower=0> [m] par2_gamma;


  // Standard quantities
  vector<lower=0>[p] sdv; // feature standard deviatons
  vector[p] means;        // feature means
}

// Parametes block (important).
parameters{
  vector[p] beta_star;           // scaled feature effect
  vector<lower=0>[m] gamma_star; // scaled BP basis effect
}

transformed parameters{
      vector[n] log_lik;              //  log-likelihood for future DIC computation
      vector[p] beta;                 // beta coefficients
      vector<lower=0>[m] gamma;       // BP gamma
      vector[n] exp_eta = exp(X * beta_star);

      // vector[p] beta_std = (beta_star - to_vector(location_beta)) ./ to_vector(scale_beta);    // reparametrized beta (due to HM dynamics)

      beta = beta_star ./ sdv;      // beta to original scale

      if(M == 2){	       // if AFT
        vector<lower=0>[n] y = time ./exp_eta;
        vector[n] y_alt = y ./ max(tau_aft);
        matrix[n, m] g2;
        matrix[n, m] G2;

        for(i in 1:n){
          for(j in 1:m){
            g2[i, j] = beta_lpdf(y_alt[i] | j, (m - j + 1));
            G2[i, j] = beta_lcdf(y_alt[i]| j, (m - j + 1));
          }
        }
        g2 = exp(g2) ./ tau_aft;
        G2 = exp(G2);

          gamma = gamma_star * exp(sum(beta_star .* means ./ sdv));
      }
      else{	// if PO or PH
          gamma = gamma_star * exp(-sum(beta_star .* means ./ sdv));
      }

      if(null == 1){ // null
          log_lik = loglik_null(gamma_star, status, X, g, G, M, rand, id, z);
      }
      else{
        if(M == 0){ // P0
            log_lik = loglik_po(exp_eta, gamma_star, status, X, g, G, rand, id, z);
        }
        else if( M == 1){ // PH
            log_lik = loglik_ph(exp_eta, gamma_star, status, X, g, G, rand, id, z);
        }
        else{ // AFT
            log_lik = loglik_aft(time, exp_eta, gamma_star, gamma, status, X, g2, G2, rand, id, z);
        }
      }
  }

// Model block (important).
model{

  if(approach == 1){ // priors
  ////// Beta Prior
      for(i in 1:p){
        if(priordist_beta[i] == 0){
          beta_star[i] ~ normal(location_beta[i], scale_beta[i]);
        }
        else{
          beta_star[i] ~ cauchy(location_beta[i], scale_beta[i]);
        }
      }
  ////// Gamma Prior
      for(i in 1:m){
        if(priordist_gamma[i] == 0){
          gamma_star[i] ~ gamma(par1_gamma[i], par2_gamma[i]);
        }
        else if(priordist[i] == 1){
          gamma_star[i] ~ inv_gamma(par1_gamma[i], par2_gamma[i]);
        }
        else{
          gamma_star[i] ~ lognormal(par1_gamma[i], par2_gamma[i]);
        }
      }
   }
     target += sum(log_lik);
 }

// Final line empty to avoid warnings.
