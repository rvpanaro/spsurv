// Functions block (optional).
// Here we implement a likelihood function.
functions{

   vector loglikph(vector beta, vector gamma, vector status, matrix Z, matrix b, matrix B, int null){

    vector[num_elements(status)] loglik;
    vector[num_elements(status)] h0;
    vector[num_elements(status)] H0;
    vector[num_elements(status)] theta;

    h0 = b*gamma;
    H0 = B*gamma;
    theta = exp(Z*beta);
    
	if(null == 1){
	loglik = log(h0) .* status - H0;
	}
	else{
	loglik = log(h0 .* theta) .* status - H0 .* theta;
	}
    return loglik;
  }

  vector loglikpo(vector beta, vector gamma, vector status, matrix Z, matrix b, matrix B, int null){

    vector[num_elements(status)] loglik;
    vector[num_elements(status)] r0;
    vector[num_elements(status)] R0;
    vector[num_elements(status)] theta;

    r0 = b * gamma;
    R0 = B * gamma;
    theta = exp(Z*beta);

	if(null == 1){
	loglik = log( (r0) ./ (1 + R0)) .* status - log(1 + R0);
	}
	else{
	loglik = log( (r0 .* theta) ./ (1 + R0 .* theta)) .* status - log(1 + R0 .* theta);
	}  
    return loglik;
  }
}

// Final line empty to avoid warnings.
