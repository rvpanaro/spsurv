survivor <- function(spbp, ...){
  UseMethod("survivor", spbp)
}

#' Survivor function calculations for Bernstein Polynomial based regression models
#'
#' @export
#' @description A method to ease survivor function computation.
#' @param n integer; sample size.
#'
#' @details sim_surv returns weibull (log-logistic) randomly
#' generated survival times. According to Collett (2003), the
#' Accelerated Failure Time emcompasses a wide variety of parametric
#' models, including weibull and log-logistic models.
#'
#' @return data.frame of `ncol(x) +2` columns in which the
#'  survival times are the response variable denoted by `y`,
#'   `status` indicates failure (0 = failure) and the features
#'   are appended to the next columns.
#'
#'@examples
#'library(spsurv)
#'
#' db <- sim_surv(100, model = "aft")
#' fitmle <- spbp(Surv(y,status) ~ x1 + x2,  model = "aft",
#'           approach = "mle", data = db)
#'
#' curve(survivor(time = x, spbp = fitmle), ylim = c(0,1), xlim = c(0,fit1$tau))
#'
#' points(x = 2, y = survivor(time = 2, fit1), pch = 19, col = 3)
#'
#' @seealso  \code{\link[spsurv]{spbp}}, \code{\link[spsurv]{sim_surv}}
#' @references
#'
#' Osman, M., & Ghosh, S. K. (2012). Nonparametric regression models for right-censored data using Bernstein polynomials. Computational Statistics & Data Analysis, 56(3), 559-573.

survivor.default <- function(time,
                             arg = list(beta = NULL, gamma = NULL),
                             newdata,
                             model = c("ph", "po", "aft"),
                             approach = c("mle", "bayes")){

  if(sum(names(arg) %in% c("beta", "gamma")) != 2)
    stop('`args` names do not match')

  ## CALL EXCEPTION HANDLING
  approach <- match.arg(approach)
  model <- match.arg(model)
  beta <- arg$beta
  gamma <- arg$gamma

  if(!is.vector(time))
    stop("time is not a vector")

  if(!is.vector(gamma))
    stop("gamma is not a vector")

  if(!is.vector(beta))
    stop("beta is not a vector")

    if(!is.data.frame(newdata))
      stop("newdata is not a data.frame")

  x <- newdata
  degree <- length(gamma)
  k <- 1:degree
  y <- time[order(time)]
  tau <- max(y)
  B <- matrix(sapply(k, function(k) pbeta(y/tau, k, degree - k + 1)), ncol = degree)
  eta <- as.matrix(x) %*% matrix(beta, ncol = 1)

  if(model == "ph"){
    H0 <- apply(B, 1, function(x){gamma %*% x})
    H <- as.vector(exp(eta)) * H0
  }
  else if(model == "po"){
    R0 <- apply(B, 1, function(x){gamma  %*% x})
    R <- as.vector(exp(eta)) * R0
    H <- -log(1 + R)
  }
  else{
    y_aft <- y / exp(eta)
    tau_aft <- max(y_aft)
    B <- matrix(sapply(k, function(k) pbeta(y_aft / tau_aft, k, degree - k + 1)), ncol = degree)
    H <- apply(B, 1, function(x){gamma %*% x})
  }
  return(exp(-H))
}

#' Survivor function calculations for an spbp object outcome
#'
#' @export
#' @description A method to call survivor.default.
#' @param n integer; sample size.
#'
#' @details sim_surv returns weibull (log-logistic) randomly
#' generated survival times. According to Collett (2003), the
#' Accelerated Failure Time emcompasses a wide variety of parametric
#' models, including weibull and log-logistic models.
#'
#' @return data.frame of `ncol(x) +2` columns in which the
#'  survival times are the response variable denoted by `y`,
#'   `status` indicates failure (0 = failure) and the features
#'   are appended to the next columns.
#'
#'@examples
#'library(spsurv)
#'
#' db <- sim_surv(100, model = "aft")
#' fit1 <- spbp(Surv(y,status) ~ x1 + x2,  model = "aft",
#'           approach = "mle", data = db)
#'
#' curve(survivor(time = x, spbp = fit1),
#' ylim = c(0,1), xlim = c(0,fit1$tau))
#'
#'points(x = 2, y = survivor.spbp(time = 2, fit1), pch = 19, col = 3)
#'
#' @seealso \code{\link[spsurv]{sim_llogis}}, \code{\link[spsurv]{spbp}}
#' @references
#'
#' Osman, M., & Ghosh, S. K. (2012). Nonparametric regression models for right-censored data using Bernstein polynomials. Computational Statistics & Data Analysis, 56(3), 559-573.

survivor.spbp <- function(spbp, newdata){
  design <- model.matrix(spbp)
  if(missing(newdata)){
    newdata <- data.frame(t(matrix(colMeans(design))))
    colnames(newdata) <- colnames(design)
  }
  if(!is.data.frame(newdata))
    stop("newdata is not a data.frame object")
  if(ncol(newdata) != ncol(design))
    stop("cols must match with `model.matrix(spbp)`")

  if(spbp$call$approach == "bayes"){
    beta <- rstan::extract(spbp$stanfit, pars = "beta_std")$beta
    gamma <- rstan::extract(spbp$stanfit, pars = "gamma_std")$gamma_std
    iter <- nrow(beta)
    ####

    s <- matrix(NA, ncol = length(spbp$y[,1]), nrow = iter)
    for(i in 1:iter){
      s[i, ] <- survivor.default(time = spbp$y[,1],
                                    arg = list(beta = beta[i,],
                                               gamma = gamma[i,]),
                                    newdata = newdata,
                                    model = spbp$call$model,
                                    approach = spbp$call$approach)
    }
  }
  else{
    beta <- spbp$coefficients[1:spbp$q]
    gamma <- spbp$coefficients[(spbp$q+1):length(spbp$coefficients)]
    ####

        s <- survivor.default(time = spbp$y[,1],
                     arg = list(beta = beta,
                                gamma = gamma),
                     newdata = newdata,
                     model = spbp$call$model,
                     approach = spbp$call$approach)
  }
  return(s)
}
