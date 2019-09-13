#' Cumulative Hazard function for Bernstein Polynomial based regression
#'
#' @export
#' @param spbp a spbp class object.
#' @param degree Bernstein Polynomial degree.
#' @param tau Real valued number greater than any time observed.
#' @param data a data.frame object.
#' @param approach Bayesian or Maximum Likelihood estimation methods, default is approach = "bayes".
#' @param model Proportional Hazards or Proportional Odds BP based regression, default is model = "ph".
#' @param priors Prior settings for the Bayesian approach.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains) or `rstan::optimizing`.
#' @return An object of class `stanfit` returned by `rstan::stan`.
#' @seealso https://mc-stan.org/users/documentation/
#' @examples
#' print(fitmle)
#' @references Osman, M. and Ghosh, S. K. (2012), “Nonparametric regression models for right-censoreddata using Bernstein polynomials,”Computational Statistics & Data Analysis, 56, 559–573.
#' @importFrom rstan stan sampling optimizing
#' @importFrom survival Surv

cumhaz <- function(t, features, coef, degree, tau = NULL, model){
  x <- features
  q <- (length(coef) - degree)
  beta <- matrix(coef[1:q], ncol = 1)
  gamma <- matrix(coef[(q + 1):length(coef)], ncol = 1)

  if(model %in% c('ph', 'po')){
    if(is.null(tau))
      tau <- max(t)

    cumhaz.f <- function(t, x, beta, gamma, degree){
      # b <- bp(time = t, m = degree, tau = tau)$b
      B <- bp(time = t, m = degree, tau = tau)$B

      if(model == 'ph'){
        return(B %*% gamma * exp(x %*% beta))
      }
      else{
        return(log(1 + B %*% exp(x %*% beta)))
      }
    }
  }
  else if( model == 'aft'){

    cumhaz.f <- function(t, x, beta, gamma, degree){
      t <- t / exp(x * beta)
      if(is.null(tau))
        tau <- max(t)

      # b <- bp(time = t, m = degree, tau = tau)$b
      B <- bp(time = t, m = degree, tau = tau)$B

      return(B %*% gamma)
    }
  }
  else{
    stop('There is no such model')
  }
  result <- sapply(t, function(x){
    apply(features, 1, cumhaz.f, t = x,
          beta = beta,
          gamma = gamma,
          degree = degree)})
  return(as.matrix(t(result)))
}

#' Survivor function for Bernstein Polynomial based regression
#'
#' @export
#' @param spbp a spbp class object.
#' @param degree Bernstein Polynomial degree.
#' @param tau Real valued number greater than any time observed.
#' @param data a data.frame object.
#' @param approach Bayesian or Maximum Likelihood estimation methods, default is approach = "bayes".
#' @param model Proportional Hazards or Proportional Odds BP based regression, default is model = "ph".
#' @param priors Prior settings for the Bayesian approach.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains) or `rstan::optimizing`.
#' @return An object of class `stanfit` returned by `rstan::stan`.
#' @seealso https://mc-stan.org/users/documentation/
#' @examples
#' print(fitmle)
#' @references Osman, M. and Ghosh, S. K. (2012), “Nonparametric regression models for right-censoreddata using Bernstein polynomials,”Computational Statistics & Data Analysis, 56, 559–573.
#' @importFrom rstan stan sampling optimizing
#' @importFrom survival Surv

surv <- function(t, features, coef, degree, tau = NULL, model){
  exp(-cumhaz(t = t, features = features, coef = coef, degree = degree,
         tau = tau, model = model))
}
