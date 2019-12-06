#' Survivor function calculations for Bernstein Polynomial based regression models
#'
#' @description A method to allow survivor function computation.
#' @param spbp
#'#' @examples
#' data("veteran") ## imports veteran dataset from survival package
#'
#' library("spsurv")
#'
#' fit <- spbp(Surv(time, status) ~ karno + factor(celltype),
#' data = veteran, approach =  "bayes", model = "po", chains = 1, iter = 1000)
#'
#' survivor(fit)
#'
#' @seealso \code{\link[spsurv]{spbp}}, \code{\link[spsurv]{sim_surv}}
#' @export survivor
#' @references
#'
#' Osman, M., & Ghosh, S. K. (2012). Nonparametric regression models for right-censored data using Bernstein polynomials. Computational Statistics & Data Analysis, 56(3), 559-573.
survivor <- function(spbp, ...) {
  UseMethod("survivor", spbp)
}

#' @return Returns the probabilities that a subject will survive beyond any given times.
#' @method survivor default
#' @export

survivor.default <- function(time,
                             arg = list(beta = NULL, gamma = NULL),
                             newdata,
                             model = c("ph", "po", "aft"),
                             approach = c("mle", "bayes")){
  if(nrow(newdata)>1){
    res <- apply(newdata, 1, survivor.calc,
                 time = time, arg = arg,
                 model = model, approach = approach)
  }
  else{
    res <- survivor.calc(time = time, arg = arg,
                         newdata = as.numeric(newdata), model = model,
                         approach = approach)
  }
  return(res)
}

#' @return Returns the probabilities that a subject will survive beyond any given times.
#' @method survivor spbp
#' @export

survivor.spbp <- function(spbp, newdata){
  if(missing(newdata)){
    design <- model.matrix(spbp)
    newdata <- data.frame(t(matrix(colMeans(design))))
    colnames(newdata) <- colnames(design)
  }
  if(!is.data.frame(newdata))
    stop("newdata is not a data.frame object")

  if(spbp$call$approach == "bayes"){

    beta <- rstan::extract(spbp$stanfit, pars = "beta_std")$beta_std
    if(ncol(newdata) != ncol(beta))
      stop("cols must match with `model.matrix(spbp)`")
    gamma <- rstan::extract(spbp$stanfit, pars = "gamma_std")$gamma_std
    iter <- nrow(beta)
    ####

    s <- matrix(NA, ncol = length(spbp$y[,1]), nrow = iter)
    for(i in 1:iter){
      s[i, ] <- survivor.default(time = spbp$y[,1],
                                    arg = list(beta = beta[i,],
                                               gamma = gamma[i,]),
                                    newdata = matrix(newdata, nrow = 1),
                                    model = spbp$call$model,
                                    approach = spbp$call$approach)
    }
  }
  else{
    beta <- spbp$coefficients[1:spbp$q]
    if(ncol(newdata) != length(beta))
      stop("cols must match with `model.matrix(spbp)`")
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

survivor.calc <- function(time,
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

  # if(!is.data.frame(newdata))
  #   stop("newdata is not a data.frame")

  x <- newdata
  degree <- length(gamma)
  k <- 1:degree
  y <- time[order(time)]
  tau <- max(y)
  B <- matrix(sapply(k, function(k) pbeta(y/tau, k, degree - k + 1)), ncol = degree)

  eta <- x %*% matrix(beta, ncol = 1)

  if(model == "ph"){
    H0 <- apply(B, 1, function(x){gamma %*% x})
    H <- as.vector(exp(eta)) * H0
  }
  else if(model == "po"){
    R0 <- apply(B, 1, function(x){gamma  %*% x})
    R <- as.vector(exp(eta)) * R0
    H <- log(1 + R)
  }
  else{
    y_aft <- as.matrix(y) / as.vector(exp(eta))
    tau_aft <- max(y_aft)
    B <- matrix(sapply(k, function(k) pbeta(y_aft / tau_aft, k, degree - k + 1)), ncol = degree)
    H <- apply(B, 1, function(x){gamma %*% x})
  }
  return(exp(-H))
}
