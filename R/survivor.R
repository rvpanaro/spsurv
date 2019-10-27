#' Survivor function calculations for Bernstein Polynomial based regression models
#'
#' @export
#' @description A method to ease survivor function computation.
#' @param spbp
#'
#' @details
#' @examples
#' data("veteran") ## imports from survival package
#' library("spsurv")
#'
#' fit <- spbp(Surv(time, status) ~ karno + factor(celltype),
#' data = veteran, approach =  "bayes", model = "po", chains = 1, iter = 1000)
#'
#' survivor(fit)
#'
#' @seealso  \code{\link[spsurv]{spbp}}, \code{\link[spsurv]{sim_surv}}
#' @references
#'
#' Osman, M., & Ghosh, S. K. (2012). Nonparametric regression models for right-censored data using Bernstein polynomials. Computational Statistics & Data Analysis, 56(3), 559-573.
survivor <- function(spbp, ...) {
  UseMethod("survivor", spbp)
}
#' @return Returns the probabilities that a subject will survive beyond any given times.
#' @rdname survivor
#' @method survivor default
#' @S3method survivor default

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

#' @return Returns the probabilities that a subject will survive beyond any given times.
#' @rdname survivor
#' @method survivor spbp
#' @S3method survivor spbp

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
