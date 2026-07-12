
functions {

  // cumulative hazard H(t|x)
  vector cumhaz(matrix Gbas, vector gamma, vector eta, int M) {
    int n = rows(Gbas);
    vector[n] B0 = Gbas * gamma;
    vector[n] H;

    if (M == 0) {                // PO
      H = log1p(B0 .* exp(eta));
    } else if (M == 1) {         // PH
      H = B0 .* exp(eta);
    } else {                     // AFT (baseline evaluated at y = log(t) - eta)
      H = B0;
    }
    return H;
  }

  // log hazard log h(t|x)
  vector log_haz(matrix gbas, matrix Gbas, vector gamma, vector eta, vector log_t, int M) {
    int n = rows(gbas);
    vector[n] b0 = gbas * gamma;
    vector[n] B0 = Gbas * gamma;
    vector[n] log_h;

    if (M == 0) {                // PO: h = b0 e^eta / (1 + B0 e^eta)
      log_h = log(b0) + eta - log1p(B0 .* exp(eta));
    } else if (M == 1) {         // PH: h = b0 e^eta
      log_h = log(b0) + eta;
    } else {                     // AFT: h(t|x) = h0(log(t) -eta) / t
      log_h = log(b0) - log_t;
    }
    return log_h;
  }

  // Pointwise log-likelihood (not stored in transformed parameters)
  vector bp_pointwise_log_lik(
    int M,
    vector status,
    vector log_time,
    matrix X,
    matrix g,
    matrix G,
    matrix P,
    vector beta,
    vector gamma
  ) {
    int n = rows(X);
    vector[n] eta = X * beta;
    vector[n] H;
    vector[n] log_h;

    if (M == 2) {
      matrix[n, cols(P)] b;
      matrix[n, cols(P)] B;
      vector[n] y = log_time - eta;
      real tau_a = min(y);
      real tau_b = max(y);
      real range = tau_b - tau_a;
      vector[n] u = (y - tau_a) / range;

      for (j in 1:cols(P)) {
        b[, j] = pow(u, j - 1);
        B[, j] = pow(u, j) / j;
      }

      b = (b * P) / range;
      B = (B * P);
      H = cumhaz(B, gamma, eta, M);

      if (min(u) < 0 || max(u) > 1) {
        log_h = rep_vector(negative_infinity(), n);
      } else {
        log_h = log_haz(b, B, gamma, eta, log_time, M);
      }
    } else {
      H = cumhaz(G, gamma, eta, M);
      log_h = log_haz(g, G, gamma, eta, log_time, M);
    }

    return -H + status .* log_h;
  }
}

// Final line empty to avoid warnings.
