#' Semiparametric survival analysis using Bernstein Polynomials
#'
#' @export
#' @param formula a Surv object with time to event, status and explanatory terms.
#' @param degree Bernstein Polynomial degree.
#' @param tau Real valued number greater than any time observed.
#' @param data a data.frame object.
#' @param approach Bayesian or Maximum Likelihood estimation methods, default is approach = "bayes".
#' @param model Proportional Hazards or Proportional Odds BP based regression, default is model = "ph".
#' @param priors Prior settings for the Bayesian approach.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains) or `rstan::optimizing`.
#' @return An object of class `stanfit` returned by `rstan::stan`.
#' @seealso https://mc-stan.org/users/documentation/
#' @examples
#' library(KMsurv)
#' data("larynx")
#'
#' #Fit a Proportional Odds BP based model to the veteran cancer data set
#' library ("survival")
#' data("veteran")
#'
#' library("spsurv")
#'
#' fitmle <− spbp(Surv(time, delta) ~ karno + factor(celltype),
#' data = veteran, approach =  "mle", model = "po")
#'
#' print(fitmle)
#' @references Osman, M. and Ghosh, S. K. (2012), “Nonparametric regression models for right-censoreddata using Bernstein polynomials,”Computational Statistics & Data Analysis, 56, 559–573.
#' @importFrom rstan stan sampling optimizing
#' @importFrom survival Surv

spbp <- function(formula, degree = NULL, tau = max(time), data,
                 approach = c("mle", "bayes"),
                 model = c("ph", "po", "aft"),
                 priors = list(shape_gamma = .01, rate_gamma = .01,
                               mean_beta = 0, sd_beta = 10),
                 hessian = TRUE, verbose = FALSE,
                 init = 0, algorithm = "LBFGS", ...) {

  ## --------------- Degree error handling ---------------
  ifelse(is.null(degree),
         {degree  <- ceiling(sqrt(nrow(data)))},
         {if(!is.integer(degree)) stop('Polynomial degree must be numeric.')}
)
  #-------------------------------------------------------------------
  model_flag <- model; approach_flag <- approach ### creates flags to save char
  model <- ifelse(match.arg(model) == "po", 0,
                  ifelse(match.arg(model) == "ph", 1, 2))
  approach <- ifelse(match.arg(approach) == "mle", 0, 1)

  ## --------------- Formula => model.frame args error handling ---------------
  # evaluate model.frame() containing the required formula
  Call <- match.call()
  aux <- match(c("formula", "data", "degree", "tau"),
               names(Call), nomatch = 0)
  if (aux[1] == 0) stop("A formula argument is required")
  if (aux[2] == 0) stop("A dataset argument is required")

  ## --------------- Approach error handling ---------------
  defaultPriors <- list(shape_gamma = .01, rate_gamma = .01, mean_beta = 0,
                        sd_beta = 10)
  ## case 1: bayes aproach w/ wrong prior spec.
  if(approach == 1 & !is.null(priors)){
    if(sum(c('shape_gamma','rate_gamma','mean_beta', 'sd_beta') %in% names(priors)) < 4) stop('Prior arguments do not match.')
  }
  ## case 2: ## case 1: bayes aproach wout/ prior spec
  else if(approach == 1 &  sum(priors %in% defaultPriors) == 4){
    message('Due to bayes approach, default priors are attributed,
            see approach in ??bpph().')
  }
  ## case 3: mle approach w/ prior spec.
  else{
    if(!is.null(priors))message('Due to mle approach priors are ignored.')
  }
  #-------------------------------------------------------------------

  ## --------------- Extra args error handling ---------------
  ##  ... arguments directly passed to `rstan::stan`, handles typos
  ## like "chans=4".

  stanArgs <- list(...)

  if (length(stanArgs)) {
    ifelse(approach == 0,
           stanformals <- c(names(formals(rstan::stan)),
                            "seed", "check_data", "sample_file",
                            ~~"algorithm", "verbose", "hessian", "as_vector",
                            "draws", "constrained", "save_iterations",
                            "refresh", "init_alpha", "tol_obj",
                            "tol_rel_obj", "tol_grad", "tol_rel_grad",
                            "tol_param", "history_size"),
           stanformals <- names(formals(rstan::stan))) #legal arg names
    aux <- pmatch(names(stanArgs), stanformals, nomatch = 0)

    if (any(aux == 0))
      stop(gettextf("Argument %s not matched", names(stanArgs)[aux==0]))
  }

  #-------------------------------------------------------------------
  temp <- Call[c(1, aux)]  # keep important args

  temp[[1L]] <- quote(stats::model.frame)  # model frame call
  temp$formula <- terms(formula, data = data);

  mf <- eval(temp, parent.frame())

  if (nrow(mf) == 0) stop("Only missing observations")
  Terms <- terms(mf)

  Y <- model.extract(mf, "response") # in general, time-to-event response
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
  type <- attr(Y, "type")
  if (type!='right' && type!='counting')
    stop(paste("Proportional hazards model doesn't support \"", type,
               "\" survival data", sep=''))
  if (length(attr(Terms, 'variables')) > 2) { # a ~1 formula has length 2
    ytemp <- terms.inner(formula)[1:2]
    xtemp <- terms.inner(formula)[-c(1,2)]
  if (any(!is.na(match(xtemp, ytemp))))
      warning("a variable appears on both the left and right sides of
              the formula")
  }
#-------------------------------------------------------------------

## --------------- Model Fitting ---------------
  ## Data
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

  features <- X
  attr(X, "assign") <- attr(model.matrix(Terms, mf), "assign")[-1]
  attr(X, "contrasts") <- attr(model.matrix(Terms, mf), "contrasts")
  xlevels <- .getXlevels(Terms, mf)
  contrasts <- attr(X, "contrasts")

  assign <- attrassign(X, Terms)
  X <-  scale(X, scale = T)
  q <- ncol(X)
  time <- as.vector(Y[,1])
  status <- as.vector(Y[,2])

  # base calculations
  base <- bp(time, m = degree, tau = tau)

  # data
  standata <- list(time = time, tau = tau, n = data.n, m = degree, q = q,
                   status = status, X = X, B = base$B, b = base$b,
                   approach = approach, M = model, null = null)

  ## Stanfit
  standata <- c(standata, priors)
  # mle
  if(approach == 0){
    stanfit <- rstan::optimizing(stanmodels$spbp, data = standata, init = init,
                                 hessian = hessian, verbose = verbose, ...)
    coef <- stanfit$par ## rescaled coefficients
    coef[1:q] <- stanfit$par[1:q] / attr(X, 'scaled:scale')
    beta <- coef[1:q]

    info <- chol2inv(-stanfit$hessian)/(attr(X, 'scaled:scale'))^2 ## rescaled fisher info

    names(beta) <- colnames(X)
    names(coef) <- c(colnames(X), paste0('gamma', 1:(length(stanfit$par)-q)))
    if(hessian == FALSE || null == 1){
      stanfit$hessian <- matrix(rep(NA, q^2), ncol = 1:q,
                                nrow = 1:q)
    }
    nulldata <- standata
    nulldata$null <- 1
    nullfit <- rstan::optimizing(stanmodels$spbp, data = nulldata, init = init,
                                 hessian = hessian, ...)
    output <- list(coefficients = coef,
                 var = info,
                 loglik = c(nullfit$value, stanfit$value),
                 linear.predictors = c(features %*% beta),
                 means = colMeans(features),
                 method = algorithm,
                 n = data.n,
                 nevent = sum(status),
                 q = q,
                 terms = Terms,
                 assign = assign,
                 wald.test = coxph.wtest(chol2inv(-stanfit$hessian), stanfit$par)$test,
                 y = Y,
                 formula = formula,
                 xlevels = xlevels,
                 contrasts = contrasts,
                 call = Call,
                 return_code = stanfit$return_code
    )
  }
  else{   # bayes
    output <- list()
    stanfit <- rstan::sampling(stanmodels$spbp, data = standata,
                               verbose = verbose, ...)
    output$stanfit <- stanfit
  }

  output$model <- model_flag
  output$approach <- approach_flag
  class(output) <- "spbp"
  return(output)
}

#' Bernstein Polynomials basis calculations
#' @export
#' @param time a vector of times to events.
#' @param m Bernstein Polynomial degree
#' @param tau must be greater than times maximum value observed, this value is used for correction purposes.
#' @param data a data.frame with variables named in the formula.
#' @return A list containing matrices b and B corresponding BP basis and corresponding tau value used to compute them.

bp <- function(time, m,  tau = NULL){
  n <- length(time)

  if( is.null(tau) ){tau <- max(time)}
  if(sum(time >= 0) != n | m < 0 | tau <  max(time) ){
    stop("tau must be greater than a any time-to-event observation.")
    }

  k <- 1:m
  b <- matrix(NA, n, m )
  B <- matrix(NA, n, m )
  y <- time/tau

  b <- sapply(k, function(k){dbeta(y, k, m - k + 1) / tau})
  B <- sapply(k, function(k) pbeta(y, k, m - k + 1) )

  # Equivalent to
  # for (i in 1:n){
  #   for(k in 1:m){
  #     b[i,k] <- dbeta(y[i], k, m - k + 1) / tau
  #     B[i,k] <- pbeta(y[i], k, m - k + 1)
  #   }
  # }
  return(list(b = b, B = B, m = m, tau = tau))
}

terms.inner <- function (x)
{
  if (class(x) == "formula")
    c(terms.inner(x[[2]]), terms.inner(x[[3]]))
  else if (class(x) == "call" && (x[[1]] != as.name("$") &&
                                  x[[1]] != as.name("["))) {
    if (x[[1]] == "+" || x[[1]] == "*" || x[[1]] == "-") {
      c(terms.inner(x[[2]]), terms.inner(x[[3]]))
    }
    else if (x[[1]] == as.name("Surv") || x[[1]] == as.name("rand"))
      unlist(lapply(x[-1], terms.inner))
    else terms.inner(x[[2]])
  }
  else (deparse(x))
}
