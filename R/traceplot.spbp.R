#' Traceplot method for fitted spbp models
#'
#' @export
#' @description Traceplot samples from a fitted \code{\link[spsurv]{spbp}} model.
#' @param spbp the result of a \code{\link[spsurv]{spbp}} fit.
#' @param ... arguments inherent from \code{\link[rstan]{traceplot}}.
#' @return see \code{\link[rstan]{traceplot}}.
#' @seealso \code{\link[spsurv]{spbp}}, \code{\link[spsurv]{stan_dens.spbp}}, \code{\link[spsurv]{extract.spbp}}
#' @examples
#' fit <- spbp(Surv(time, status) ~ age + sex, approach = "bayes", data = lung)
#' traceplot(fit)
#' @importFrom rstan traceplot

traceplot.spbp <-
  function(spbp, ...){
    if(spbp$call$approach == "bayes")
      rstan:::traceplot(object = spbp$stanfit, ...)
    else{
      "not applicable, change approach to 'bayes' to traceplot MCMC chains"
    }
  }

