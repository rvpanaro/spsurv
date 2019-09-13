#' Gradient for Bernstein Polynomial based regression
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

grad <- function(t, features, coef, degree, tau = NULL, model){

  x <- features
  q <- (length(coef) - degree)
  beta <- matrix(coef[1:q], ncol = 1)
  gamma <- matrix(coef[(q + 1):length(coef)], ncol = 1)

  if(model == 'ph'){
    if(is.null(tau))
      tau <- max(t)

    dbetaj <- function(j, t,  x, beta, gamma, degree){
      B <- bp(time = t, m = degree, tau = tau)$B
      -(B %*% gamma) * (x[-j] %*% matrix(beta[-j, ], ncol = 1) + beta[j]) * exp((x %*% beta)-(B %*% gamma))
    }

    dgammak <- function(k, t, x, beta, gamma, degree){
      B <- bp(time = t, m = degree, tau = tau)$B
      -(B[-k] %*% matrix(gamma[-k, ], ncol = 1) + B[k]) * exp((x %*% beta)-(B %*% gamma))
    }
  }
  else if(model == 'po'){
    if(is.null(tau))
      tau <- max(t)

    dbetaj <- function(j, t,  x, beta, gamma, degree){
      B <- bp(time = t, m = degree, tau = tau)$B
      -(1/(1 + (B %*% gamma) * exp(x %*% beta))^2) * (x[-j] %*% matrix(beta[-j, ], ncol = 1) + beta[j]) * exp(x %*% beta)
    }

    dgammak <- function(k, t, x, beta, gamma, degree){
      B <- bp(time = t, m = degree, tau = tau)$B
      -(1/(1 + (B %*% gamma) * exp(x %*% beta))^2) * -(B[-k] %*% matrix(gamma[-k, ], ncol = 1) + B[k]) * exp(x %*% beta)
    }
  }
  else{
    t <- t / exp(x * beta)
    if(is.null(tau))
      tau <- max(t)

    dbetaj <- function(j, t,  x, beta, gamma, degree){
      b <- bp(time = t, m = degree, tau = tau)$b
      B <- bp(time = t, m = degree, tau = tau)$B
      (b %*% gamma) * t *(x[-j] %*% matrix(beta[-j, ], ncol = 1) + beta[j]) * exp(-(x %*% beta)-(B %*% gamma))
    }
    dgammak <- function(k, t, x, beta, gamma, degree){
      b <- bp(time = t, m = degree, tau = tau)$b
      B <- bp(time = t, m = degree, tau = tau)$B
      -(B[-k] %*% matrix(gamma[-k, ], ncol = 1) + B[k]) * exp(-(B %*% gamma))
      }
  }
  grad = list()
  for(j in 1:q){
    grad[[j]] <- sapply(t, function(x){
      apply(features, 1, dbetaj, t = x,
            j = j,
            beta = beta,
            gamma = gamma,
            degree = degree)})
  }
  for(k in 1:degree){
    grad[[(q + k)]] <- sapply(t, function(x){
      apply(features, 1, dgammak, t = x,
            k = k,
            beta = beta,
            gamma = gamma,
            degree = degree)})
  }

  result <- list()
  if(ncol(features)>1){
    result <- cbind(sapply(grad,'['))
  }
  else{
    for(i in 1:ncol(features)){
      for(j in 1:length(coef)){
        result[[i]][, j] <- sapply(grad,'[[', j)
      }
    }
  }
 return(result)
}
