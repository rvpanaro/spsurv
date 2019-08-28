// Functions block (optional).
// Here we implement a likelihood function.
functions{

   vector loglikph(vector beta, vector gamma, vector status, matrix X, matrix b, matrix B, int null){

    vector[num_elements(status)] loglik;
    vector[num_elements(status)] h0;
    vector[num_elements(status)] H0;
    vector[num_elements(status)] theta;

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

  vector loglikpo(vector beta, vector gamma, vector status, matrix X, matrix b, matrix B, int null){

    vector[num_elements(status)] loglik;
    vector[num_elements(status)] r0;
    vector[num_elements(status)] R0;
    vector[num_elements(status)] theta;

    r0 = b * gamma;
    R0 = B * gamma;
    theta = exp(X*beta);

	if(null == 1){
	loglik = log( (r0) ./ (1 + R0)) .* status - log(1 + R0);
	}
	else{
	loglik = log( (r0 .* theta) ./ (1 + R0 .* theta)) .* status - log(1 + R0 .* theta);
	}
    return loglik;
 }

 // matrix bp(vector time, int m, real tau, vector phi){
 //
 //     //variable declaration
 //     int n = num_elements(time);
 //     vector[n] y = time ./ phi;
 //     real tau_alt = max(tau ./ phi);
 //     matrix [n, m] basis[2];
 //     matrix [n, m] b;
 //     matrix [n, m] B;
 //
 //
 //     for (i in 1:n){
 //       for(k in 1:m){
 //         b[i, k] = exp(beta_lpdf(y[i] | k, m - k + 1) / tau_alt);
 //         B[i, k] = beta_cdf(y[i] , k, m - k + 1);
 //      }
 //      basis[1] = b;
 //      basis[2] = B;
 //
 //     return basis;
 //    }
 //   }

 vector loglikaft(vector time, vector beta, vector gamma, vector status, matrix X,
    matrix b, matrix B, real tau, int null){

    vector[num_elements(status)] loglik;
    vector[num_elements(status)] h0;
    vector[num_elements(status)] H0;
    vector[num_elements(status)] phi;
    matrix[num_elements(status), num_elements(gamma)] basis[2];
    basis[1] = b;
    basis[2] = B;
    h0 = basis[1]*gamma;
    H0 = basis[2]*gamma;
    phi = exp(X * beta);

  	if(null == 1){
  	  loglik = log(h0) .* status - H0;
	  }
	  else{
  	  // basis = bp(time, num_elements(gamma), tau, phi);
// 	    h0 = basis[1] * gamma;
//       H0 = basis[2] * gamma;
		  loglik = log(h0 ./ phi) .* status - H0;
	  }
    return loglik;
 }
}
// Final line empty to avoid warnings.
