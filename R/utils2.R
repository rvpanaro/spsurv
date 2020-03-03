#---------------------------------------------
#' Generic S3 method model.matrix
#' @aliases model.matrix
#' @export
#' @param spbp a fitted model object
#' @param ... further arguments passed to or from other methods.
#' @return the matrix of explanatory variables of a data set.
#'
model.matrix <- function(spbp, ...) UseMethod("model.matrix")

#' Model.matrix method for fitted spbp models
#'
#' @aliases model.matrix.spbp
#' @rdname model.matrix-methods
#' @method model.matrix spbp
#' @export
#' @export model.matrix
#' @description Model.matrix of a fitted \code{\link[spsurv]{spbp}} model.
#' @param object an object of class `spbp`, see \code{\link[spsurv]{spbp}}.
#' @param ... arguments inherent from \code{\link[stats]{model.matrix}}.
#' @return The explanatory variables matrix.
#' @seealso \code{\link[spsurv]{spbp}}, \code{\link[stats]{model.matrix}}
#' @examples
#' fit <- spbp(Surv(time, status) ~ age + sex, approach = "bayes", data = lung)
#' model.matrix(fit)

model.matrix.spbp <-
  function(spbp, data = eval(spbp$call$data, envir = parent.frame())){
    model.matrix(as.formula(spbp$call$formula), data = data)[, -1]
  }


#---------------------------------------------
#' Generic S3 method traceplot
#' @aliases traceplot
#' @export
#' @param spbp a fitted model object
#' @param ... further arguments passed to or from other methods.
#' @return the traceplot of a MCMC chain.
#'
traceplot <- function(spbp, ...) UseMethod("traceplot")

#' Traceplot method for fitted spbp models
#'
#' @aliases traceplot.spbp
#' @rdname traceplot-methods
#' @method traceplot spbp
#' @export
#' @export traceplot
#' @description Traceplot of a Bayesian fit \code{\link[spsurv]{spbp}}.
#' @param spbp an object of class `spbp` result of a \code{\link[spsurv]{spbp}} fit.
#' @param ... arguments inherent from \code{\link[rstan]{traceplot}}.
#' @return see \code{\link[rstan]{traceplot}}.
#' @seealso \code{\link[spsurv]{spbp}}, \code{\link[spsurv]{stan_dens.spbp}}, \code{\link[spsurv]{extract.spbp}}
#' @examples
#' fit <- spbp(Surv(time, status) ~ age + sex, approach = "bayes", data = lung)
#' traceplot(fit)
#' @importFrom rstan traceplot

traceplot.spbp <-
  function(spbp, pars = c("beta", "gamma"), ...){
    if(spbp$call$approach == "bayes")
      rstan:::traceplot(object = spbp$stanfit, pars = pars, ...)
    else{
      warning("not applicable, change approach to 'bayes' to get MCMC chains density plots")
    }
  }

#---------------------------------------------
#' Generic S3 method extract
#' @aliases extract
#' @export
#' @param spbp a fitted model object
#' @param ... further arguments passed to or from other methods.
#' @return extract the MCMC chain values of a Bayesian fit.
#'
extract <- function(spbp, ...) UseMethod("extract")

#' Extract method for fitted spbp models
#'
#' @aliases extract.spbp
#' @rdname extract-methods
#' @method extract spbp
#' @export
#' @export extract
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
  function(spbp, pars = c("beta", "gamma"), ...){
    if(spbp$call$approach == "bayes")
      rstan:::extract(object = spbp$stanfit, pars = pars, ...)
    else{
      warning("not applicable, change approach to 'bayes' to get MCMC chains density plots")
    }
  }

#---------------------------------------------
#' Generic S3 method extract
#' @aliases stan_dens
#' @export
#' @param spbp a fitted model object
#' @param ... further arguments passed to or from other methods.
#' @return the density plot of a MCMC chain.
#'
stan_dens <- function(spbp, ...) UseMethod("stan_dens")

#' Density plots method for fitted spbp models
#'
#' @aliases stan_dens.spbp
#' @rdname stan_dens-methods
#' @method stan_dens spbp
#' @export
#' @export stan_dens
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
  function(spbp, pars = c("beta", "gamma"), ...){
    if(spbp$call$approach == "bayes")
      rstan:::stan_dens(object = spbp$stanfit, pars = pars, ...)
    else{
      warning("not applicable, change approach to 'bayes' to get MCMC chains density plots")
    }
  }


