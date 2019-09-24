residuals.spbp <- function(object, type=c("martingale", "deviance")){

  model_flag <- object$model
  approach_flag <- object$approach
  ni <- object$linear_predictors
  gamma <- object$coef

  if(model_flag == "ph" && approach_flag == "mle"){
    H_est =

    if(type == "martingale"){}
    else {}


  }
  else if(model_flag == "po" && approach_flag == "mle"){

  }
  else if(model_flag == "aft" && approach_flag == "mle"){

  }
  else{
    stop('No residuals for bayesian models yet.')
  }

  type <- if.else(match.arg(type), 0, 1)
}
