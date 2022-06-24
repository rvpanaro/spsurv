// Functions block (optional).
// Here we implement a likelihood function and other useful functions.
functions{

 vector loglik_null(vector gamma_star, vector status,
   matrix g, matrix G, int n, int M, int rand, int [] id, vector z){

    vector[n] log_lik;
    vector[n] bp;
    vector[n] BP;

    bp = g * gamma_star;
    BP = G * gamma_star;

    if(rand == 1){ //gamma
      bp = bp .* z[id];
      BP = BP .* z[id];
    }
    else if( rand == 2) { //additive
      bp = bp .* exp(z[id]);
      BP = BP .* exp(z[id]);
    }

    if(M == 0){
        log_lik = (log(bp) - log(1 + BP)) .* status - log(1 + BP); // prop odds
    }
    else{
        log_lik = log(bp) .* status - BP; // aft or prop hazards
    }
    return log_lik;
  }

 vector loglik_ph(vector exp_eta, vector gamma_star, vector status,
 matrix g, matrix G, int n, int rand, int [] id, vector z){

      vector[n] log_lik;
      vector[n] h0;
      vector[n] H0;

      h0 = g * gamma_star;
      H0 = G * gamma_star;

    if(rand == 1){ //gamma
      h0 = h0 .* z[id];
      H0 = H0 .* z[id];
    }
    else if( rand == 2) { //additive
      h0 = h0 .* exp(z[id]);
      H0 = H0 .* exp(z[id]);
    }

    log_lik = (log(h0) + log(exp_eta)) .* status - (H0 .* exp_eta);
 return log_lik;
}

 vector loglik_po(vector exp_eta, vector gamma_star, vector status,
  matrix g, matrix G, int n, int rand, int [] id, vector z){

      vector[n] log_lik;
      vector[n] r0;
      vector[n] R0;

      r0 = g * gamma_star;
      R0 = G * gamma_star;

    if(rand == 0){
      r0 = r0 .* z[id];
      R0 = R0 .* z[id];
    }
    else if( rand == 2) { //additive
      r0 = r0 .* exp(z[id]);
      R0 = R0 .* exp(z[id]);
    }

	  log_lik = (log((r0 .* exp_eta)) - log(1 + R0 .* exp_eta)) .* status - log(1 + R0 .* exp_eta);
 return log_lik;
 }

 vector loglik_aft(vector time,  vector exp_eta, vector gamma_star, vector gamma, vector status,
    matrix g, matrix G, int n, int rand, int [] id, vector z){

      vector[n] log_lik;
      vector[n] h0;
      vector[n] H0;

      h0 = g * gamma_star;
      H0 = G * gamma;

    if(rand == 1){ //gamma multiplicative
      h0 = g * gamma_star .* z[id];
      H0 = G * gamma_star .* z[id];
    }
    else if( rand == 2) { //additive
      h0 = g * gamma_star .* exp(z[id]);
      H0 = G * gamma_star .* exp(z[id]);
    }

	  log_lik =  ((log(h0)-log(exp_eta)) .* status) - H0;
 return log_lik;
 }
}
// Final line empty to avoid warnings.
