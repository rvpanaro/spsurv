#' Bernstein Polynomials based Proportional Hazards
#'
#' @export
#' @description Fits Bernstein Polynomial based Proportional Odds model to lung cancer data.
#' @param formula a Surv object with time to event, status and explanatory terms.
#' @param degree Bernstein Polynomial degree.
#' @param tau Real valued number greater than any time observed.
#' @param data a data.frame object.
#' @param approach Bayesian or Maximum Likelihood estimation methods, default is approach = "bayes".
#' @param model Proportional Hazards or Proportional Odds BP based regression, default is model = "ph".
#' @param priors Prior settings for the Bayesian approach.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains) or `rstan::optimizing`.
#' @return An object of class `stanfit` returned by `rstan::stan`.
#' @seealso \url{https://mc-stan.org/users/documentation/}
#' @examples
#'
#' data("veteran") ## imports from survival package
#' library("spsurv")
#'
#' fit <- spbp(Surv(time, status) ~ karno + factor(celltype),
#' data = veteran, approach =  "bayes", model = "po", chains = 1, iter = 1000)
#'
#' print(fit)
#' summary(fit)
#' @references Osman, M. and Ghosh, S. K. (2012), “Nonparametric regression models for right-censoreddata using Bernstein polynomials,”Computational Statistics & Data Analysis, 56, 559–573.
#' @importFrom rstan stan sampling optimizing
#' @importFrom survival Surv frailty


bpph <- function(formula, data, ...){
  spbp.default(formula = formula,
               data = data,
               model = "ph",
               approach = "mle",
               init = init,
               ...)
}

#' Bernstein Polynomials based Proportional Odds
#'
#' @export
#' @description Fits Bernstein Polynomial based Proportional Odds model to lung cancer data.
#' @param formula a Surv object with time to event, status and explanatory terms.
#' @param degree Bernstein Polynomial degree.
#' @param tau Real valued number greater than any time observed.
#' @param data a data.frame object.
#' @param approach Bayesian or Maximum Likelihood estimation methods, default is approach = "bayes".
#' @param model Proportional Hazards or Proportional Odds BP based regression, default is model = "ph".
#' @param priors Prior settings for the Bayesian approach.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains) or `rstan::optimizing`.
#' @return An object of class `stanfit` returned by `rstan::stan`.
#' @seealso \url{https://mc-stan.org/users/documentation/}
#' @examples
#'
#' data("veteran") ## imports from survival package
#' library("spsurv")
#'
#' fit <- spbp(Surv(time, status) ~ karno + factor(celltype),
#' data = veteran, approach =  "bayes", model = "po", chains = 1, iter = 1000)
#'
#' print(fit)
#' summary(fit)
#' @references Osman, M. and Ghosh, S. K. (2012), “Nonparametric regression models for right-censoreddata using Bernstein polynomials,”Computational Statistics & Data Analysis, 56, 559–573.
#' @importFrom rstan stan sampling optimizing
#' @importFrom survival Surv frailty

bppo <- function(formula, data, init = 0,...){
  spbp.default(formula = formula,
               data = data,
               model = "po",
               approach = "mle",
               init = init,
               ...)
}

#' Bernstein Polynomials based Accelerated Failure Time
#'
#' @export
#' @description Fits Bernstein Polynomial based Proportional Odds model to lung cancer data.
#' @param formula a Surv object with time to event, status and explanatory terms.
#' @param degree Bernstein Polynomial degree.
#' @param tau Real valued number greater than any time observed.
#' @param data a data.frame object.
#' @param approach Bayesian or Maximum Likelihood estimation methods, default is approach = "bayes".
#' @param model Proportional Hazards or Proportional Odds BP based regression, default is model = "ph".
#' @param priors Prior settings for the Bayesian approach.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains) or `rstan::optimizing`.
#' @return An object of class `stanfit` returned by `rstan::stan`.
#' @seealso \url{https://mc-stan.org/users/documentation/}
#' @examples
#'
#' data("veteran") ## imports from survival package
#' library("spsurv")
#'
#' fit <- spbp(Surv(time, status) ~ karno + factor(celltype),
#' data = veteran, approach =  "bayes", model = "po", chains = 1, iter = 1000)
#'
#' print(fit)
#' summary(fit)
#' @references Osman, M. and Ghosh, S. K. (2012), “Nonparametric regression models for right-censoreddata using Bernstein polynomials,”Computational Statistics & Data Analysis, 56, 559–573.
#' @importFrom rstan stan sampling optimizing
#' @importFrom survival Surv frailty

bpaft <- function(formula, data, ...){
  spbp.default(formula = formula,
               data = data,
               model = "aft",
               approach = "mle",
               init = init,
               ...)
}
