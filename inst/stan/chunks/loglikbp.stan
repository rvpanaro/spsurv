// Functions block (optional).
// Here we implement a likelihood function.
functions{
  real new_beta(real x, real alpha, real beta){
    real prob;
    prob = inc_beta(alpha, beta, 1)*(x^(alpha-1))*((1-x)^(beta-1));
    return prob;
  }

   vector loglikph(vector beta, vector gamma, vector status, matrix X, matrix b,
    matrix B, int null){
      int n = num_elements(status);
      vector[n] loglik;
      vector[n] h0;
      vector[n] H0;
      vector[n] theta;

      h0 = b*gamma;
      H0 = B*gamma;
      theta = exp(X*beta);

  	if(null == 1){
	    loglik = log(h0) .* status - H0;
	  }
	  else{
	   loglik = log(h0 .* theta) .* status - H0 .* theta;
	  }
    return loglik;
 }

  vector loglikpo(vector beta, vector gamma, vector status, matrix X, matrix b,
    matrix B, int null){
      int n = num_elements(status);
      vector[n] loglik;
      vector[n] r0;
      vector[n] R0;
      vector[n] theta;

    r0 = b * gamma;
    R0 = B * gamma;
    theta = exp(X * beta);

	  if(null == 1){
	    loglik = log( (r0) ./ (1 + R0)) .* status - log(1 + R0);
	  }
	  else{
	    loglik = log( (r0 .* theta) ./ (1 + R0 .* theta)) .* status - log(1 + R0 .* theta);
	  }
    return loglik;
 }

 vector loglikaft(vector time, vector beta, vector gamma, vector status, matrix X,
    matrix b, matrix B, int null){
    // declaration
    int n = num_elements(status);
    int m = num_elements(gamma);
    vector[n] loglik;
    vector[n] h0;
    vector[n] H0;
    vector[n] phi = exp(X * beta);
    vector[n] y = time ./ phi;
    real tau_alt = max(time ./ phi);
    matrix[n, m] b_alt;
    matrix[n, m] B_alt;


  	if(null == 1){
  	  h0 = b * gamma;
  	  H0 = B * gamma;
  	  loglik = log(h0) .* status - H0;
	  }
	  else{
      for (i in 1:n){
        for(k in 1:m){
          b_alt[i, k] = new_beta((y[i]/ tau_alt), k, (m - k + 1)) / tau_alt;
          B_alt[i, k] = beta_cdf((y[i]/ tau_alt), k, (m - k + 1));
        }
      }
      h0 = b_alt * gamma;
      H0 = B_alt * gamma;
	    loglik = log(h0 ./ phi) .* status - H0;
	  }
    return loglik;
 }
}
// Final line empty to avoid warnings.
