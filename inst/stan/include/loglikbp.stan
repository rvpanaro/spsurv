// Functions block (optional).
// Here we implement a likelihood function.
functions{

 vector loglik_null(vector beta, vector nu, vector status, matrix X, matrix b,
  matrix B, int M, int dist, int[] id, vector z){
    int n = num_elements(status);
    vector[n] loglik;
    vector[n] bp;
    vector[n] BP;
    vector[n] theta = exp(X * beta);

    if(dist == 0){
      bp = b * exp(nu);
      BP = B * exp(nu);
    }
    else{
      bp = b * exp(nu + z[id]);
      BP = B * exp(nu + z[id]);
    }
    if(M == 0){
        loglik = (log(bp) - log(1 + BP)) .* status - log(1 + BP); // prop odds
    }
    else{
        loglik = log(bp) .* status - BP; // aft or prop hazards
    }
    return loglik;
  }

 vector loglik_ph(vector beta, vector nu, vector status, matrix X, matrix b,
  matrix B, int dist, int[] id, vector z){
    int n = num_elements(status);
    vector[n] loglik;
    vector[n] h0;
    vector[n] H0;
    vector[n] theta = exp(X * beta);

    if(dist == 0){
      h0 = b * exp(nu);
      H0 = B * exp(nu);
    }
    else{
      h0 = b * exp(nu + z[id]);
      H0 = B * exp(nu + z[id]);
    }
    loglik = (log(h0) + log(theta)) .* status - (H0 .* theta);
 return loglik;
}

 vector loglik_po(vector beta, vector nu, vector status, matrix X, matrix b,
  matrix B, int dist, int[] id, vector z){
    int n = num_elements(status);
    vector[n] loglik;
    vector[n] r0;
    vector[n] R0;
    vector[n] theta = exp(X * beta);

    if(dist == 0){
      r0 = b * exp(nu);
      R0 = B * exp(nu);
    }
    else{
      r0 = b * exp(nu + z[id]);
      R0 = B * exp(nu + z[id]);
    }
	  loglik = (log((r0 .* theta)) - log(1 + R0 .* theta)) .* status - log(1 + R0 .* theta);
 return loglik;
 }

 vector loglik_aft(vector time, vector beta, vector nu, vector status, matrix X,
    matrix b, matrix B, int dist, int[] id, vector z){
    // declaration
  int n = num_elements(status);
  int m = num_elements(exp(nu));
  vector[n] loglik;
  vector[n] h0;
  vector[n] H0;
  vector[n] phi = exp(X * beta);
  vector[n] y = time ./ phi;
  real tau = max(y) + 0.000001; // prevents z becoming 1
  vector[n] y_alt = y ./ tau;
  matrix[n, m] b2;
  matrix[n, m] B2;

	 for (i in 1:n){
    for(k in 1:m){
     b2[i, k] = exp(beta_lpdf(y_alt[i]| k, (m - k + 1))) / tau;
     B2[i, k] = beta_cdf(y_alt[i], k, (m - k + 1));
    }
   }
    if(dist == 0){
      h0 = b2 * exp(nu);
      H0 = B2 * exp(nu);
    }
    else{
      h0 = b2 * exp(nu + z[id]);
      H0 = B2 * exp(nu + z[id]);
    }
	  loglik =  ((log(h0)-log(phi)) .* status) - (H0);
 return loglik;
 }
}
// Final line empty to avoid warnings.
