
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
  vector log_haz(matrix gbas, matrix Gbas, vector gamma, vector eta, vector log_t,  int M) {
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
}

// Final line empty to avoid warnings.
