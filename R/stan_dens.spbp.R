stan_dens.spbp <-
  function(spbp, ...){
    if(spbp$call$approach == "bayes")
      rstan:::stan_dens(object = spbp$stanfit, ...)
    else{
      "not applicable, change approach to 'bayes' to get MCMC chains density plots"
    }
  }

