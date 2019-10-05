traceplot.spbp <-
  function(spbp, ...){
    if(spbp$call$approach == "bayes")
      rstan:::traceplot(object = spbp$stanfit, ...)
    else{
      "not applicable, change approach to 'bayes' to traceplot MCMC chains"
    }
  }

