#' Density plots method for fitted spbp models
#'
#' @export
#' @method stan_dens spbp
#' @description Posterior density of samples from a fitted \code{\link[spsurv]{spbp}} model.
#' @param spbp the result of a \code{\link[spsurv]{spbp}} fit.
#' @param ... arguments inherent from \code{\link[rstan]{stan_dens}}.
#' @return see \code{\link[rstan]{stan_dens}}.
#' @seealso \code{\link[spsurv]{spbp}}, \code{\link[spsurv]{traceplot.spbp}}, \code{\link[spsurv]{extract.spbp}}
#' @examples
#' fit <- spbp(Surv(time, status) ~ age + sex, approach = "bayes", data = lung)
#' stan_dens(fit)
#' @importFrom rstan stan_dens
#'
stan_dens.spbp <-
  function(spbp, ...){
    if(spbp$call$approach == "bayes")
      rstan:::stan_dens(object = spbp$stanfit, ...)
    else{
      "not applicable, change approach to 'bayes' to get MCMC chains density plots"
    }
  }

