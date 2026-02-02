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
  
}

parameters {
  vector[p] beta;                  // beta 1-SD unit (needs to be transformed to the original scale)
  vector<lower=0>[m] gamma;        // baseline BP coefficients on normalized design matrix
}

transformed parameters {
  vector[n] eta   = X * beta;      // = sum_j beta_j (x_ij - mu_j)
  real alpha      = 0.0;

  vector[n] H;
  vector[n] log_h;

  if (M == 2) {
    matrix[n,m] b;
    matrix[n,m] B;

    vector[n] y = log_time - eta;

    real tau_a = min(log_time) - max(eta);// virtual 0 in exp scale
    real tau_b = max(log_time) - min(eta);                      

    real range = (tau_b - tau_a);
    vector[n] u = fmin(fmax((y - tau_a) / fmax(range, 1e-12), 1e-12), 1 - 1e-12);

    for (j in 1:m) {
      b[,j] = pow(u, j - 1);       // n-vector
      B[,j] = pow(u, j) / j;       // integral in u
    }

    b = (b * P) / range;
    B = (B * P);
    
    H     = cumhaz(B, gamma, eta, M);
    log_h = log_haz(b, B, gamma, eta, log_time, M);
    
  } else {
    
    alpha = sum(beta .* means ./ sdv);
    H     = cumhaz(G, gamma, eta, M);     // No matrix copies here: just use g, G directly
    log_h = log_haz(g, G, gamma, eta, log_time, M);
    
  }
  
  vector[n] log_lik = -H + status .* log_h;
}

model {
  target += sum(log_lik);          // likelihood

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

// Final line empty to avoid warnings.
