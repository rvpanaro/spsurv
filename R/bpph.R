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
                 a_gamma = .01, b_gamma = .01,
                 m_beta = 0, S_beta = 100,
                 chains = 1, ...) {

  Call <- match.call()
  approach = ifelse(match.arg(approach) == "frequentist", 0, 1)

  #Step 1
  ## We want to pass any ... args to `rstan::sampling`, but not pass things
  ##  like "chans=4" where someone just made a typo.  The use of ...
  ##  is simply to allow things like "iter=5000" with easier typing
  extraArgs <- list(...)
  if (length(extraArgs)) {
    stanargs <- names(formals(rstan::stan)) #legal arg names
    indx <- pmatch(names(extraArgs), stanargs, nomatch=0L)
    if (any(indx==0L))
      stop(gettextf("Argument %s not matched",
                    names(extraArgs)[indx==0L]), domain = NA)
  }

  # Step 2
  # create a call to model.frame() that contains the formula (required)
  #  and any other of the relevant optional arguments
  # then evaluate it in the proper frame
  indx <- match(c("formula", "m", "tau", "optimize",
                  "a_gamma", "b_gamma", "m_beta", "S_beta", "data"),
                names(Call), nomatch=0)
  if (indx[1] == 0) stop("A formula argument is required")
  temp <- Call[c(1,indx)]  # only keep the arguments we wanted
  temp[[1L]] <- quote(stats::model.frame)  # change the function called

  temp$formula <- if(missing(data)) terms(formula)
  else              terms(formula, data=data);

  mf <- eval(temp, parent.frame());

  if (nrow(mf) ==0) stop("No (non-missing) observations")
  Terms <- terms(mf)

  Y <- model.extract(mf, "response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
  type <- attr(Y, "type")
  if (type!='right' && type!='counting')
    stop(paste("Cox model doesn't support \"", type,
               "\" survival data", sep=''))
  data.n <- nrow(Y)   #remember this before any time transforms
  if (length(attr(Terms, 'variables')) > 2) { # a ~1 formula has length 2
    ytemp <- terms.inner(formula)[1:2]
    xtemp <- terms.inner(formula)[-c(1,2)]
    if (any(!is.na(match(xtemp, ytemp))))
      warning("a variable appears on both the left and right sides of the formula")
  }
  ## ends Survival package

  Z = model.matrix(Terms, mf)
  time <- as.matrix(Y)[,1]
  status <- as.matrix(Y)[,2]

  base <- bp(time, m = m, tau = tau)

  standata <- list(n = data.n, m = m, q = ncol(Z[,-1]),
                   status = as.vector(status), Z = Z[,-1], B = base$B, b = base$b,
                   approach = approach)
  # frequentist
  if(approach == 0){
    standata$a_gamma = a_gamma
    standata$b_gamma = b_gamma
    standata$m_beta = m_beta
    standata$S_beta = S_beta

          stanfit <- rstan::optimizing(stanmodels$bpph, data = standata,...)

    summary <- summary(stanfit)$summary
    return(print(summary))
  }
  # bayesian
  else{
    stanfit <- rstan::sampling(stanmodels$bpph, data = standata, chains = chains,  ...)
    return(stanfit)
  }
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
