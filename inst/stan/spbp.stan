// include Functions block.
#include /include/loglikbp.stan

data {
  int<lower=1> n;
  int<lower=1> m;
  int<lower=1> p;

  int<lower=0,upper=1> approach;   // 0 = freq, 1 = Bayes
  int<lower=0,upper=2> M;          // 0 = PO; 1 = PH; 2 = AFT

  vector<lower=0,upper=1>[n] status;
  vector<lower=0>[n] time;

  matrix[n,p] X;
  matrix[n,m] g;
  matrix[n,m] G;
  matrix[m,m] P;                   // Power basis

  int<lower=0> priordist_beta[p];
  vector[p] location_beta;
  vector<lower=0>[p] scale_beta;

  int<lower=0> priordist_gamma[m];
  vector<lower=0>[m] location_gamma;
  vector<lower=0>[m] scale_gamma;

  vector[p] means;
  vector<lower=0>[p] sdv;

}

transformed data {

  vector[n] log_time = log(time);  // log-time scale
  // Minimum feasible divisor to avoid 0/0 when range = 0; ~1.5 * machine epsilon (double)
  real min_range = 1e-12;

}

parameters {
  vector[p] beta;                  // beta 1-SD unit (needs to be transformed to the original scale)
  vector<lower=0>[m] gamma;        // baseline BP coefficients on normalized design matrix
}

transformed parameters {
  real alpha = (M == 2) ? 0.0 : sum(beta .* means ./ sdv);
}

model {
  target += sum(
    bp_pointwise_log_lik(M, status, log_time, X, g, G, P, beta, gamma)
  );

  if (approach == 1) {             // priors
    for (j in 1:p) {
      if (priordist_beta[j] == 0)
        target += normal_lpdf(beta | location_beta[j], scale_beta[j]) + log(sdv[j]);
      else
        target += logistic_lpdf(beta | location_beta[j], scale_beta[j]) + log(sdv[j]);
    }

    for (j in 1:m) {
      if (priordist_gamma[j] == 0)
        gamma[j] ~ lognormal(location_gamma[j], scale_gamma[j]);
      else {
        gamma[j] ~ normal(location_gamma[j], scale_gamma[j]);
      }
    }
  }
}

generated quantities {
  vector[approach == 1 ? n : 0] log_lik;
  if (approach == 1) {
    log_lik = bp_pointwise_log_lik(
      M, status, log_time, X, g, G, P, beta, gamma
    );
  }
}

// Final line empty to avoid warnings.
