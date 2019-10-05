extract.spbp <-
  function(spbp, ...){
    if(spbp$call$approach == "bayes")
      rstan:::extract(object = spbp$stanfit, ...)
    else{
      "not applicable, change approach to 'bayes' to extract MCMC chains"
    }
  }
