spbp <- function(formula, ...){
  UseMethod("spbp", formula)
}

#' Semiparametric survival analysis using Bernstein Polynomials
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

spbp.default <-
  function(formula, degree, data,
            approach = c("mle", "bayes"),
            model = c("ph", "po", "aft"),
            priors = list(beta = c("normal(0,4)"),
                         gamma = "lognormal(0,4)"),
           scale = TRUE,
           ...){

  # ---------------Definitions + error handling  ---------------
  ## tau degree

  if(missing(degree)){
    degree <- ceiling(sqrt(nrow(data)))
  }

  ## Call
  Call <- match.call();

  ## model
  if(length(model) == 3){model_flag = "ph";}
  else{model_flag <- model;}

  model <- ifelse(match.arg(model) == "po", 0,
                  ifelse(match.arg(model) == "ph", 1, 2))

  ## approach
  approach_flag <- match.arg(approach) ### saves string input
  approach <- ifelse(approach_flag == "mle", 0, 1)

  handler1() ## error handling #1

  ## terms
  temp <- Call[c(1, aux)] # keep important args
  temp[[1L]] <- quote(stats::model.frame) # model frame call
  special <- c("frailty", "frailty.gamma", "frailty.gaussian", "frailty.t")
  temp$formula <- terms(formula, special, data = data)
  temp$formula <- terms(formula, data = data);

  ## frailty (id, distribution, column)
  handler2()

  ## Priors
  suppressMessages(handler3())

  ## stanArgs
  stanArgs <- list(...)
  handler4()

  ## Model Frame

  mf <- eval(temp, parent.frame())
  Terms <- terms(mf)
  Y <- model.extract(mf, "response") # time-to-event response
  type <- attr(Y, "type")
  handler5()

  # ---------------  Data declaration + definitions ---------------

  ## + sample size + labels
  data.n <- nrow(Y)
  labels <- attributes(temp$formula)$term.labels
  null <- 0

  if (length(labels) > 1){
    X <-  model.matrix(Terms, mf)[, -1]
  }
  else if(length(labels) == 1){
    X <- as.matrix(model.matrix(Terms, mf)[, -1], ncol = data.n)
    colnames(X) <- labels
  }
  else{
    X <- as.matrix(rep(0, data.n), ncol = data.n)
    colnames(X) <- "non-parametric"
    null <- 1
  }

  ## time + status + features
  features <- X
  attr(X, "assign") <- attr(model.matrix(Terms, mf), "assign")[-1]
  attr(X, "contrasts") <- attr(model.matrix(Terms, mf), "contrasts")
  xlevels <- .getXlevels(Terms, mf)
  contrasts <- attr(X, "contrasts")

  assign <- attrassign(X, Terms)

  q <- ncol(X)
  time <- as.vector(Y[,1])
  tau <- max(time)
  status <- as.vector(Y[,2])

  if(scale == T){
    X <- scale(X)
    std <- attr(X, "scaled:scale")

    ## rescaled coefficients (correction)
    if(model == 2){ ## AFT
      wsum <- exp(sum(attr(X, "scaled:center")/ attr(X, "scaled:scale")))
    }
    else{
      wsum <- exp(-sum(attr(X, "scaled:center")/ attr(X, "scaled:scale")))
    }
  }
  else{
    std <- rep(1, q)
    wsum <- 1
  }

  ## base calculations
  base <- bp.basis(time, degree = degree, tau = tau)

  ## priors to num

  priordist <- sapply(priordist,
                      function(x){
                        switch(x,
                               "normal" = 0,
                               "gamma" = 1,
                               "inv_gamma" = 2,
                               "lognormal" = 3)})

  priordist_beta <- sapply(priordist_beta,
                          function(x){switch(x,
                          "normal" = 0,
                          "cauchy" = 1)})
  ## Recycling

  priordist_beta <- array(priordist_beta, dim = q)
  location_beta <- array(location_beta, dim = q)
  scale_beta <- array(scale_beta, dim = q)

  ## standata
  standata <- list(time = time,
                   tau = tau,
                   n = data.n,
                   m = base$degree,
                   q = q,
                   status = status,
                   X = X,
                   B = base$B,
                   b = base$b,
                   approach = approach,
                   M = model,
                   null = null,
                   id = rep(1, data.n),
                   dist = dist,
                   z = rep(0, data.n),
                   priordist = as.numeric(priordist),
                   priorpars = as.numeric(priorpars),
                   priordist_beta = as.numeric(priordist_beta),
                   location_beta = as.numeric(location_beta),
                   scale_beta = as.numeric(scale_beta),
                   std  = std,
                   wsum = wsum
  )

  # --------------- Fit  ---------------

  output <- list()
  if(approach == 0){
    spbp.mle(standata = standata, ...)
  }
  else{
    spbp.bayes(standata = standata, ...)
  }
}

spbp.mle <-
  function(standata,
           init = 0,
           hessian = TRUE,
           verbose = FALSE,
           ...){

  e <- parent.frame()
  approach_flag <- get("approach_flag", envir = e)
  model_flag <- get("model_flag", envir = e)
  features <- get("features", envir = e)
  X <- get("X", envir = e)
  Y <- get("Y", envir = e)
  q <- get("q", envir = e)
  degree <- get("degree", envir = e)
  frailty_idx <- get("frailty_idx", envir = e)
  null <- get("null", envir = e)
  data.n <- get("data.n", envir = e)
  status <- get("status", envir = e)
  Terms <- get("Terms", envir = e)
  xlevels <- get("xlevels", envir = e)
  tau <- get("tau", envir = e)
  Call <- get("Call", envir = e)

  if(!is.null(frailty_idx)){
    standata$X <- X[, -frailty_idx]
    message("Frailty ignored, change approach to `bayes` for frailty estimation.")
  }

  stanfit <-
    rstan::optimizing(stanmodels$spbp,
                             data = standata,
                             init = init,
                             hessian = hessian,
                             verbose = verbose,
                             ...)
  len <- length(stanfit$par)
  ## stanfit coefficients (beta, nu)
  coef <- stanfit$par[(len-q - degree+1):len]
  ## regression estimates
  beta <- coef[1:q]

  ## rescaled hessian matrix
  hess <- stanfit$hessian
  hess[1:q, 1:q] <- stanfit$hessian[1:q, 1:q]

  names(beta) <- colnames(X)
  names(coef) <- c(names(beta),
                   paste0("gamma", 1:(degree))
  )

  ## singular matrices handler
  # if(det(-hess) == 0)
  #   stop("Optimizing hesssian matrix is singular!")

  ## rescaled fisher info
  info <- blockSolve(hess, q)
  diag(info) <- abs(diag(info))

  if(hessian == FALSE || null == 1){
    stanfit$hessian <- matrix(rep(NA, q^2),
                              ncol = 1:q,
                              nrow = 1:q)
  }
  nulldata <- standata
  nulldata$null <- 1
  nullfit <- rstan::optimizing(stanmodels$spbp,
                               data = nulldata,
                               init = init,
                               hessian = hessian,
                               ...)

  output <- list(coefficients = coef,
                 var = info,
                 loglik = c(nullfit$value, stanfit$value),
                 linear.predictors = c(features %*% beta),
                 means = colMeans(features),
                 method = "optimizing",
                 n = data.n,
                 nevent = sum(status),
                 q = q,
                 terms = Terms,
                 assign = assign,
                 wald.test = coxph.wtest(info[1:q, 1:q], beta)$test,
                 y = Y,
                 formula = formula,
                 xlevels = xlevels,
                 contrasts = contrasts,
                 return_code = stanfit$return_code,
                 tau = tau,
                 call = Call)
  output$call$approach <- approach_flag
  output$call$model <- model_flag

  class(output) <- "spbp"
  message('Priors are ignored due to mle approach .')

  return(output)
}

spbp.bayes <- function(standata,
                       hessian = TRUE,
                       verbose = FALSE,
                       chains = 1,
                       ...){
  e <- parent.frame()
  approach_flag <- get("approach_flag", envir = e)
  model_flag <- get("model_flag", envir = e)
  Call <- get("Call", envir = e)
  features <- get("features", envir = e)
  q <- get("q", envir = e)
  degree <- get("degree", envir = e)
  frailty_idx <- get("frailty_idx", envir = e)
  X <- get("X", envir = e)
  Y <- get("Y", envir = e)
  dist <- get("dist", envir = e)

  # bayes
  output <- list(y = Y)

  if(dist == 0){
    output$stanfit <- rstan::sampling(stanmodels$spbp,
                                      data = standata,
                                      verbose = verbose,
                                      chains = chains,
                                      ...)

    output$pmode <- rstan::optimizing(stanmodels$spbp,
                                      data = standata,
                                      verbose = verbose)$par[1:(q+degree)]
  }
  else{
    standata$X <- X[, -frailty_idx]

    output$stanfit <- rstan::sampling(stanmodels$spbp_frailty,
                                      data = standata,
                                      verbose = verbose,
                                      chains = chains,
                                      ...)

    output$pmode <- rstan::optimizing(stanmodels$spbp_frailty,
                                      data = standata,
                                      verbose = verbose)$par[1:(q+degree)]
  }
  output$loo <- loo::loo(loo::extract_log_lik(output$stanfit))
  output$waic <- loo::waic(loo::extract_log_lik(output$stanfit))
  output$call <- Call
  output$call$approach <- approach_flag
  output$call$model <- model_flag
  class(output) <- "spbp"
  return(output)
}
