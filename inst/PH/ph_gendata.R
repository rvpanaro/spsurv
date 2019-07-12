exp_surv_gen <- function(n, beta, design_matrix, p = .8){
  X <- cbind(1, design_matrix);

  if( length(beta) != (dim(X)[2] ))
    return(warning("Beta length doesn't match covariates column dimension!"))

  u <- runif(n = n, 0, 1)
  failt <- - log(u) / exp(X %*% beta)

  ## Censorship percentage p= P[T<C]
  lambda_t <- exp(X %*% beta)
  lambda_c <- ((1-p) * lambda_t)/p  # p == lambda_t/(lambda_c + lambda_t)

  c <- rexp(n = n, rate = lambda_c)

  ## Censorship indicator
  delta <- ifelse( failt < c, 1, 0)

  ## t right censored
  t <- ifelse( failt < c, failt, c)

  data= data.frame(time = t, status = delta, X[,-1])
  return(data)
}
