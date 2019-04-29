#' Semiparametric survival analysis using Bernstein Polynomials
#'
#' @export
#' @param formula a formula object, with time-to-event data on the left side of ~ and explanatory terms on the right.
#' @param degree Bernstein Polynomial degree.
#' @param tau must be greater than times maximum value observed, used for correction purposes.
#' @param approach bayes or mle methods, default is approach = "bayes".
#' @param model proportional hazards or porportional odds base, default is model = ph
#' @param data a data.frame with variables named in the formula.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains), see https://mc-stan.org/users/documentation/.
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'

spbp <- function(formula, degree = NULL, tau = NULL, data,
                 approach = c("bayes", "mle"),
                 model = c("ph", "po"),
                 priors = list(shape_gamma = .01, rate_gamma = .01, mean_beta = 0, sd_beta = 10),
                 verbose = FALSE, ...) {

  ## --------------- Degree error handling ---------------
  ifelse(is.null(degree),
         {degree  <- ceiling(sqrt(nrow(data)))},
         {if(!is.integer(degree)) stop('Polynomial degree must be numeric.')})
  #-------------------------------------------------------------------

  model = ifelse(match.arg(model) == "po", 0, 1)

  ## --------------- Approach error handling ---------------
  approach = ifelse(match.arg(approach) == "mle", 0, 1)

  ## case 1: bayes aproach w/ wrong prior spec.
  if(approach == 1 & !is.null(priors)){
    if(sum(c('shape_gamma','rate_gamma','mean_beta', 'sd_beta') %in% names(priors)) < 4) stop('Prior arguments do not match.')
  }
  ## case 2: ## case 1: bayes aproach wout/ prior spec
  else if(approach == 1 & is.null(priors)){
    warning('Due to bayes approach, default priors are attributed, see approach in ??bpph().')
  }
  ## case 3: mle approach w/ prior spec.
  else{
    if(!is.null(priors))warning('Due to mle approach priors are ignored.')
  }
  #-------------------------------------------------------------------

  ## --------------- Extra args error handling ---------------
  ##  ... arguments directly passed to `rstan::stan`, handles typos like "chans=4".
  stanArgs <- list(...)

  if (length(stanArgs)) {
    stanformals <- names(formals(rstan::stan)) #legal arg names
    aux <- pmatch(names(stanArgs), stanformals, nomatch = 0)

    if (any(aux == 0))
      stop(gettextf("Argument %s not matched", names(stanArgs)[aux==0]))
  }
  #-------------------------------------------------------------------

  ## --------------- Formula => model.frame args error handling ---------------
  # evaluate model.frame() containing the required formula
  Call <- match.call()

  aux <- match(c("formula", "data", "degree", "tau"), names(Call), nomatch = 0)

  if (aux[1] == 0) stop("A formula argument is required")
  if (aux[2] == 0) stop("A dataset argument is required")
  temp <- Call[c(1, aux)]  # keep important args

  temp[[1L]] <- quote(stats::model.frame)  # model frame call
  temp$formula <- terms(formula, data = data);
  mf <- eval(temp, parent.frame());

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
      warning("a variable appears on both the left and right sides of the formula")
  }
#-------------------------------------------------------------------

## --------------- Model Fitting ---------------
  ## Data
  data.n <- nrow(Y)
  if (length(attr(Terms, 'variables')) > 2){
    Z = model.matrix(Terms, mf)[,-1]
  }
  else{
    Z = as.matrix(rep(0, data.n), nrow = data.n)
  }

  Z = scale(Z, scale = F)
  time <- as.vector(Y[,1])
  status <- as.vector(Y[,2])

  base <- bp(time, m = degree, tau = tau)

  standata <- list(n = data.n, m = degree, q = ncol(Z),
                   status = status, Z = Z, B = base$B, b = base$b,
                   approach = approach, M = model)

  ## Stanfit
  standata <- do.call(c, list(standata, priors))
  # mle
  if(approach == 0){
    stanfit <- rstan::optimizing(stanmodels$spbp, data = standata, init = 0, ...)
  }
  else{   # bayes
    stanfit <- rstan::sampling(stanmodels$spbp, data = standata, ...)
  }
  return(stanfit)
}

#' Bernstein Polynomials basis calculations
#'
#' @export
#' @param time a vector of times to events.
#' @param m Bernstein Polynomial degree
#' @param tau must be greater than times maximum value observed, this value is used for correction purposes.
#' @param data a data.frame with variables named in the formula.
#' @return A list containing matrices b and B corresponding BP basis and corresponding tau value used to compute them.
#'

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

#' A function taken from the survival library - it is not exported from there hence a local copy
#'
#'@param x an R object
#'@return a character vector
#'@importFrom survival Surv
#'@keywords internal

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
