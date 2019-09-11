// Functions block (optional).
// Here we implement a likelihood function.
functions{

 vector loglik_null(vector beta, vector gamma, vector status, matrix X, matrix b,
  matrix B, int M){
    int n = num_elements(status);
    vector[n] loglik;
    vector[n] bp = b * gamma;
    vector[n] BP = B * gamma;
    vector[n] theta = exp(X * beta);

    if(M == 0){
        loglik = (log(bp) - log(1 + BP)) .* status - log(1 + BP); // prop odds
    }
    else{
        loglik = log(bp) .* status - BP; // aft or prop hazards
    }
    return loglik;
  }

 vector loglik_ph(vector beta, vector gamma, vector status, matrix X, matrix b,
  matrix B){
    int n = num_elements(status);
    vector[n] loglik;
    vector[n] h0= b * gamma;
    vector[n] H0 = B * gamma;
    vector[n] theta = exp(X * beta);

    loglik = (log(h0) + log(theta)) .* status - (H0 .* theta);
 return loglik;
}

 vector loglik_po(vector beta, vector gamma, vector status, matrix X, matrix b,
  matrix B){
    int n = num_elements(status);
    vector[n] loglik;
    vector[n] r0 = b * gamma;
    vector[n] R0 = B * gamma;
    vector[n] theta = exp(X * beta);

	  loglik = (log((r0 .* theta)) - log(1 + R0 .* theta)) .* status - log(1 + R0 .* theta);
 return loglik;
 }

 vector loglik_aft(vector time, vector beta, vector gamma, vector status, matrix X,
    matrix b, matrix B){
    // declaration
  int n = num_elements(status);
  int m = num_elements(gamma);
  vector[n] loglik;
  vector[n] h0;
  vector[n] H0;
  vector[n] phi = exp(X * beta);
  vector[n] y = time ./ phi;
  real tau = max(y) + 0.000001; // prevents z becoming 1
  vector[n] z = y ./ tau;
  matrix[n, m] b2;
  matrix[n, m] B2;

	 for (i in 1:n){
    for(k in 1:m){
     b2[i, k] = exp(beta_lpdf(z[i]| k, (m - k + 1))) / tau;
     B2[i, k] = beta_cdf(z[i], k, (m - k + 1));
    }
   }
    h0 = b2 * gamma;
    H0 = B2 * gamma;
	  loglik =  ((log(h0)-log(phi)) .* status) - (H0);
 return loglik;
 }
}
// Final line empty to avoid warnings.
