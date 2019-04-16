// Functions block (optional).
// Here we implement a likelihood function.
functions{

   vector loglikbp(vector beta, vector gamma, vector status, matrix Z, matrix b, matrix B, int M ){

    vector[num_elements(status)] loglik;
    vector[num_elements(status)] bp;
    vector[num_elements(status)] BP;
    vector[num_elements(status)] eta;

    bp = b * gamma;
    BP = B * gamma;
    eta = Z * beta;


    if(M == 1){ // proportional hazards
      loglik = log(bp .* exp(eta)) .* status - BP .* exp(eta);
    }
    else{ // proportional odds
      loglik = (eta + log(bp)  - log(exp(eta) + BP)) .* status +
              log(1 + exp(eta) .* BP);
    }
    return loglik;
  }
}

// Final line empty to avoid warnings.
