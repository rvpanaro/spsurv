#' Extract method for fitted spbp models
#'
#' @export
#' @method extract spbp
#' @description Extract samples from a fitted \code{\link[spsurv]{spbp}} model.
#' @param spbp an object of class `spbp` result of a \code{\link[spsurv]{spbp}} fit.
#' @param ... arguments inherent from \code{\link[rstan]{extract}}.
#' @return see \code{\link[rstan]{extract}}.
#' @seealso \code{\link[spsurv]{spbp}}, \code{\link[spsurv]{stan_dens.spbp}}, \code{\link[spsurv]{traceplot.spbp}}
#' @examples
#' fit <- spbp(Surv(time, status) ~ age + sex, approach = "bayes", data = lung)
#' extract(fit)
#' @importFrom rstan extract

extract.spbp <-
  function(spbp, ...){
    if(spbp$call$approach == "bayes")
      rstan:::extract(object = spbp$stanfit, ...)
    else{
      "not applicable, change approach to 'bayes' to extract MCMC chains"
    }
  }
