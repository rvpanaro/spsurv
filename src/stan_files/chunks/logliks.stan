// Functions block (optional).
// Here we implement a likelihood function.
functions{
   vector loglikbp(vector beta, vector gamma, vector status, matrix Z, matrix b, matrix B){

    vector[num_elements(status)] loglik;
    vector[num_elements(status)] h;	
    vector[num_elements(status)] H;
    vector[num_elements(status)] eta;
    h = b*gamma;
    H = B*gamma;
    eta = Z*beta;	
	
    loglik = log(h .* exp(eta)) .* status - H .* exp(eta);
    return loglik;
  }
}

// Final line empty to avoid warnings.
