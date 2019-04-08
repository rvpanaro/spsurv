#' Bernstein Polynomials base risk estimation in Proportional hazards model
#'
#' @export
#' @param formula a formula object, with time-to-event data on the left side of ~ and explanatory terms on the right.
#' @param m Bernstein Polynomial degree
#' @param tau must be greater than times maximum value observed, this value is used for correction purposes.
#' @param optimize Bayesian or frequentist inference method, by default optimize = FALSE.
#' @param data a data.frame with variables named in the formula.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'

bpph <- function(formula, m = ceiling(sqrt(nrow(data))), tau = NA, data,
                 approach = c("bayesian", "frequentist"),
                 prior = NULL,
                 chains = 1,
                 verbose = FALSE, ...) {

  ## Approach error handling
  approach = ifelse(match.arg(approach) == "frequentist", 0, 1)

  if(approach == 1 & !is.null(prior)){
    if(sum(c('a_gamma','b_gamma','m_beta', 'S_beta') %in% names(prior)) < 4) stop('Prior arguments do not match!')
    else  prior <- as.list(prior)
  }
  else if(approach == 1 & is.null(prior)){
      prior <- list(a_gamma = .01, b_gamma = .01, m_beta = 0, S_beta = 100)
      warning('Due to bayesian approach, non informative priors are attributed, see approach in ??bpph().')
  }
  else{
    if(!is.null(prior))warning('Due to frequentist approach, inputed priors are ignored.')
    prior <- list(a_gamma = 0, b_gamma = 0, m_beta = 0, S_beta = 0)
  }

#-------------------------------------------------------------------
  ##  ... arguments directly passed to `rstan::stan`, handles typos like "chans=4".
  stanArgs <- list(...)
  stanArgs <- list(chains = 1, iter = 1000)

  if (length(stanArgs)) {
    stanformals <- names(formals(rstan::stan)) #legal arg names
    aux <- pmatch(names(stanArgs), stanformals, nomatch = 0)

    if (any(aux == 0))
      stop(gettextf("Argument %s not matched", names(stanArgs)[aux==0]))
  }

  Call <- match.call()

  # evaluate model.frame() containing the required formula
  aux <- match(c("formula", "data", "m", "tau"), names(Call), nomatch = 0)

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
##---------------------------------------------------------- ends Survival package

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

  base <- bp(time, m = m, tau = tau)

  standata <- list(n = data.n, m = m, q = ncol(Z),
                   status = status, Z = Z, B = base$B, b = base$b,
                   approach = approach)

  ## Stan options

  # avoid recompilations
  rstan::rstan_options(auto_write = TRUE)

  # run different chains in parallel.
  options(mc.cores = parallel::detectCores())

  ## Stanfit
  pars = c("gamma", "beta")
  print(standata)

  model = rstan::stan_model('src/stan_files/bpph.stan')
  standata <- do.call(c, list(standata, prior))
  print(standata)

  # frequentist
  if(approach == 0){

    stanfit <- rstan::optimizing(model, data = standata, ...)
  }
  # bayesian
  else{
      stanfit <- rstan::sampling(model, data = standata, chains = chains,
                               pars = pars, ...)
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

bp <- function(time, m,  tau = NA){
  n <- length(time)

  if( is.na(tau) ){tau <- max(time)}
  if(sum(time >= 0) != n | m < 0 | tau <  max(time) ){
    stop("tau is less than a time to event observation, it must be greater!")
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
