survivor <- function(spbp, ...){
  UseMethod("survivor", spbp)
}

#' Survivor function calculations for an spbp object outcome
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

survivor.spbp <- function(spbp, time = spbp$y[, 1]){

  y <- time
  coef <- spbp$coefficients
  beta <- matrix(coef[1:spbp$q], ncol = 1)
  log_gamma <- coef[(spbp$q + 1):length(coef)]
  degree <- length(log_gamma);   k <- 1:degree
  tau <- spbp$tau
  x <- colMeans(eval(spbp$call$data)[, -c(1,2)])

  B <- matrix(sapply(k, function(k) pbeta(y/tau, k, degree - k + 1)), ncol = degree)
  H0 <- apply(B, 1, function(x){exp(log_gamma) %*% x})
  eta <- x %*% beta
  H <- as.vector(exp(eta)) * H0
  return(exp(-H)[order(y)])
}
