// Functions block (optional).
// Here we implement a likelihood function.
functions{

   vector loglikph(vector beta, vector gamma, vector status, matrix Z, matrix b, matrix B){

    vector[num_elements(status)] loglik;
    vector[num_elements(status)] h0;
    vector[num_elements(status)] H0;
    vector[num_elements(status)] eta;

    h0 = b*gamma;
    H0 = B*gamma;
    eta = Z*beta;

    loglik = log(h0 .* exp(eta)) .* status - H0 .* exp(eta);
    return loglik;
  }

  // vector loglikpo(vector beta, vector gamma, vector status, matrix Z, matrix b, matrix B){
  //
  //   vector[num_elements(status)] loglik;
  //   // vector[num_elements(status)] r0;
  //   // vector[num_elements(status)] R0;
  //   vector[num_elements(status)] h0;
  //   vector[num_elements(status)] H0;
  //   vector[num_elements(status)] eta;
  //
  //   r0 = b*gamma;
  //   R0 = B*gamma;
  //
  //   // h0 = b*gamma;
  //   // H0 = B*gamma;
  //   eta = Z*beta;
  //
  //   loglik = (eta + log(r0)  - log(exp(eta) + R0)) .* status +
  //   log(1 + exp(eta) .* R0);
  //   // loglik = log(h0 .* exp(eta)) .* status - H0 .* exp(eta);
  //   return loglik;
  // }

  vector loglikpo(vector beta, vector gamma, vector status, matrix Z, matrix b, matrix B){

    vector[num_elements(status)] loglik;
    vector[num_elements(status)] r0;
    vector[num_elements(status)] R0;
    vector[num_elements(status)] eta;

    r0 = b*gamma;
    R0 = B*gamma;
    eta = Z*beta;

    loglik = (eta + log(r0)  - log(exp(eta) + R0)) .* status +
              log(1 + exp(eta) .* R0);
    return loglik;
  }
}

// Final line empty to avoid warnings.