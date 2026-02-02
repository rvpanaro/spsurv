#' Semiparametric Survival Analysis Using Bernstein Polynomial
#'
#' Fits Bernstein Polynomial based Proportional regression to survival data.
#'
#' @title spbp: The BP Based Survival Analysis Function
#' @param formula a Surv object with time to event, status and explanatory terms.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains) or `rstan::optimizing`.
#' @seealso \code{\link[spsurv]{spbp.default}}
#' @examples
#'
#' library("spsurv")
#' data("veteran") ## imports from survival package
#'
#' fit_mle <- spbp(Surv(time, status) ~ karno + factor(celltype),
#'   data = veteran, model = "po"
#' )
#' summary(fit_mle)
#'
#' fit_bayes <- spbp(Surv(time, status) ~ karno + factor(celltype),
#'   data = veteran, model = "po", approach = "bayes",
#'   cores = 1, iter = 300, chains = 1,
#'   priors = list(
#'     beta = c("normal(0,5)"),
#'     gamma = "halfnormal(0,5)"
#'   )
#' )
#'
#' summary(fit_bayes)
#' @rdname spbp
#' @export spbp
#' @importFrom magrittr `%>%`
#' @seealso  \code{\link[spsurv]{spbp.default}},  \code{\link[spsurv]{bpph}},  \code{\link[spsurv]{bppo}}, \code{\link[spsurv]{bpaft}}, \url{https://mc-stan.org/users/documentation/}
#' @return An object of class 'spbp'.

spbp <- function(formula, ...) {
  UseMethod("spbp", formula)
}

#' @title spbp: The BP Based Semiparametric Survival Analysis Function
#' @param formula a Surv object with time to event, status and explanatory terms
#' @param degree Bernstein Polynomial degree
#' @param data a data.frame object
#' @param approach Bayesian or Maximum Likelihood estimation methods, default is approach = "bayes"
#' @param model Proportional Hazards or Proportional Odds BP based regression, default is model = "ph"
#' @param priors prior settings for the Bayesian approach; `normal` or `cauchy` for beta; `lognormal` or `loglogistic` for gamma (BP coefficients)
#' @param scale logical; indicates whether to center and scale the data
#' @param ... further arguments passed to or from other methods
#' @param cores number of core threads to use
#' @return An object of class \code{spbp}
#' @method spbp default
#' @export
#' @importFrom rstan stan sampling optimizing
#' @importFrom survival Surv frailty
#' @importFrom MASS ginv
#' @importFrom mnormt pd.solve
#' @importFrom loo waic loo
#' @importFrom coda  HPDinterval
#' @importFrom stats .getXlevels as.formula contrasts dbeta density dist formula median model.extract pbeta pchisq printCoefmat qnorm rlogis rnorm rweibull sd terms
#'

spbp.default <-
  function(formula,
           degree,
           data,
           approach = c("mle", "bayes"),
           model = c("ph", "po", "aft"),
           priors = list(
             beta = c("normal(0,4)"),
             gamma = c("lognormal(0,4)"),
             frailty = c("gamma(0.01,0.01)")
           ),
           cores = min(parallel::detectCores() - 1, 4),
           scale = TRUE,
           verbose = FALSE, chains = 4,
           ...) {
    # ---------------Definitions + error handling  ---------------
    Call <- match.call()

    if (length(model) == 3) {
      model_flag <- "ph"
    } else {
      model_flag <- model
    }

    model <- ifelse(match.arg(model) == "po", 0,
      ifelse(match.arg(model) == "ph", 1, 2)
    )

    approach_flag <- match.arg(approach)
    approach <- ifelse(approach_flag == "mle", 0, 1)

    handler1() ## error-handling nº1 -- BP degree handling

    temp <- Call[c(1, aux)] # keep important args
    temp[[1L]] <- quote(stats::model.frame) # model frame call
    special <- c("frailty", "frailty.gamma", "frailty.gaussian", "frailty.t")
    temp$formula <- terms(formula, special, data = data)
    temp$formula <- terms(formula, data = data)
    stanArgs <- list(...)

    mf <- eval(temp, parent.frame())
    Terms <- terms(mf)
    Y <- model.extract(mf, "response")
    type <- attr(Y, "type")

    time <- as.vector(Y[, 1])
    tau <- max(time)
    status <- as.vector(Y[, 2])

    if (missing(degree)) degree <- ceiling(nrow(data)^(0.5))

    handler2() ## error-handling nº2 -- Frailty (id, distribution, column)
    handler3() ## error-handling nº3 -- Priors
    handler4() ## error-handling nº5 -- Model Frame

    # ---------------  Data declaration + definitions ---------------
    data.n <- nrow(Y) ## + sample size + labels
    labels <- attributes(temp$formula)$term.labels
    null <- 0

    if (length(labels) > 0) {
      X <- model.matrix(Terms, mf)[, -1]

      if (is.null(ncol(X))) {
        X <- as.matrix(X)
        colnames(X) <- labels
      }
    } else {
      X <- as.matrix(rep(0, data.n), ncol = 1)
      colnames(X) <- "non-parametric"
      null <- 1
    }

    features <- X
    attr(X, "assign") <- attr(model.matrix(Terms, mf), "assign")[-1]
    attr(X, "contrasts") <- attr(model.matrix(Terms, mf), "contrasts")
    xlevels <- .getXlevels(Terms, mf)
    contrasts <- attr(X, "contrasts")

    assign <- attrassign(X, Terms)
    p <- ncol(X)

    if (scale & (null == 0)) {
      X <- scale(X, center = TRUE)
      means <- array(attr(X, "scaled:center"))
      sdv <- array(attr(X, "scaled:scale"))
    } else {
      means <- array(0, p)
      sdv <- array(1, p)
    }

    if (!(degree %% 1 == 0)) {
      stop("Polynomial degree must be integer.")
    }

    bp <- bp.basis(time, degree = degree, tau = tau)
    pw <- pw.basis(degree = degree)

    priordist_beta <- sapply(priordist_beta, function(x) {
      switch(x,
        "normal" = 0,
        "logistic" = 1
      )
    })

    priordist_gamma <- sapply(priordist_gamma, function(x) {
      switch(x,
        "lognormal" = 0,
        "halfnormal" = 1
      )
    })

    priordist_frailty <- sapply(priordist_frailty, function(x) {
      switch(x,
        "gamma" = 0,
        "invgamma" = 1,
        "lognormal" = 2
      )
    })

    priordist_beta <- array(priordist_beta, dim = p)
    location_beta <- array(as.numeric(location_beta), dim = p)
    scale_beta <- array(as.numeric(scale_beta), dim = p)

    priordist_gamma <- array(priordist_gamma, dim = degree)
    location_gamma <- array(as.numeric(location_gamma), dim = degree)
    scale_gamma <- array(as.numeric(scale_gamma), dim = degree)

    priordist_frailty <- array(priordist_frailty, dim = 1)
    par1_frailty <- array(as.numeric(par1_frailty), dim = 1)
    par2_frailty <- array(as.numeric(par2_frailty), dim = 1)

    standata <- list(
      n = data.n,
      m = bp$degree,
      p = p,
      tau = tau,
      approach = approach,
      rand = rand,
      M = model,
      status = status,
      id = rep(1, data.n),
      z = rep(0, data.n),
      time = time,
      X = X,
      g = bp$g,
      G = bp$G,
      P = pw,
      priordist_beta = priordist_beta,
      location_beta = location_beta,
      scale_beta = scale_beta,
      priordist_gamma = priordist_gamma,
      location_gamma = location_gamma,
      scale_gamma = scale_gamma,
      priordist_frailty = priordist_frailty,
      par1_frailty = par1_frailty,
      par2_frailty = par2_frailty,
      means = means,
      sdv = sdv
    )

    # --------------- Fit  ---------------
    if (approach == 0) {
      tryCatch(
        expr = spbp.mle(standata = standata, hessian = T, verbose = verbose, ...),
        error = function(e) {
          warning(e)
          return(NaN)
        }
      )
    } else {
      tryCatch(
        expr = spbp.bayes(standata = standata, hessian = F, verbose = verbose, chains = chains, ...),
        error = function(e) {
          warning(e)
          return(NaN)
        }
      )
    }
  }

#'
#' @export
spbp.mle <-
  function(standata, hessian = TRUE, verbose = FALSE, ...) {
    e <- parent.frame() # "sourcing" the parent.frame
    vnames <- objects(, envir = e) # variable names in parent frame

    for (n in vnames) assign(n, get(n, e))

    if (!is.null(frailty_idx)) {
      standata$X <- X[, -frailty_idx]
      message("Frailty ignored, change approach to `bayes` for frailty estimation.")
    }

    c <- 0
    stanfit <- list(return_code = 70)

    while (stanfit$return_code != 0 && c < 15) {
      # Run optimizing() with args from the list
      stanfit <- suppressWarnings(
        do.call(rstan::optimizing, c(list(
          object  = stanmodels$spbp,
          data    = standata,
          hessian = hessian,
          verbose = verbose
        ), ...))
      )
      c <- c + 1
    }

    if (stanfit$return_code != 0) {
      warning("In .local(object, ...) : non-zero return code in optimizing")
    }

    nms <- names(stanfit$par)
    alpha <- stanfit$par[nms == "alpha"] ## intercept
    coef <- stanfit$par[startsWith(nms, "beta[")] / sdv ## regression estimates
    gamma <- stanfit$par[startsWith(nms, "gamma[")] / exp(alpha) ## bernstein coefficients

    names(coef) <- colnames(X)
    names(gamma) <- paste0("gamma[", 1:degree, "]")

    if (hessian == FALSE) {
      stanfit$hessian <- matrix(rep(NA, p^2), ncol = 1:p, nrow = 1:p)
    }

    nulldata <- standata
    nulldata$X <- matrix(0, ncol = ncol(X), nrow = data.n)

    nullfit <- tryCatch(rstan::optimizing(stanmodels$spbp, data = nulldata, hessian = hessian, ...),
      error = function(e) {
        return(list(value = NULL))
      }
    )

    output <- list(
      coefficients = coef,
      bp.param = gamma,
      alpha = alpha,
      hessian = stanfit$hessian,
      loglik = c(nullfit$value, stanfit$value),
      features = features,
      n = data.n,
      nevent = sum(status),
      terms = Terms,
      y = Y,
      formula = formula,
      xlevels = xlevels,
      contrasts = contrasts,
      return_code = stanfit$return_code,
      call = Call,
      means = means,
      sdv = sdv,
      standata = standata,
      tau_a = 0,
      tau_b = tau,
      stanfit = stanfit
    )

    if (model_flag == "aft") {
      output$tau_a <- min(log(Y[, 1])) - max(features %*% coef)
      output$tau_b <- max(log(Y[, 1])) - min(features %*% coef)
    }

    if (null) {
      output$coefficients <- NULL
      output$hessian <- stanfit$hessian[-1, -1]
      output$features <- NULL
    }

    output$call$approach <- approach_flag
    output$call$model <- model_flag

    class(output) <- c("spbp")
    message("Priors are ignored because the MLE approach is used.")

    return(output)
  }

#'
#' @export
spbp.bayes <- function(standata, hessian = FALSE, verbose = FALSE, chains = 1, ...) {
  e <- parent.frame() # "sourcing" the parent.frame
  vnames <- objects(, envir = e) # variable names in parent frame

  for (n in vnames) assign(n, get(n, e))

  if (rand > 0) standata$X <- X[, -frailty_idx]

  stanfit <-
    rstan::sampling(stanmodels$spbp, data = standata, verbose = verbose, chains = chains, cores = cores, ...)

  message("\nExtracting posterior draws (this may take a moment)...\n")

  alpha_vec <- rstan::extract(
    stanfit,
    permuted = TRUE,
    pars = "alpha",
    inc_warmup = FALSE
  )$alpha

  beta_mat <- rstan::extract(
    stanfit,
    permuted = TRUE,
    pars = "beta",
    inc_warmup = FALSE
  )$beta

  gamma_mat <- rstan::extract(stanfit, pars = "gamma")$gamma

  posterior <-
    list(
      alpha = alpha_vec,
      beta = t(apply(beta_mat, 1, function(x) x / sdv)),
      gamma = apply(gamma_mat, 2, function(x) x / exp(alpha_vec)),
      log_lik = rstan::extract(
        stanfit,
        permuted = TRUE,
        pars = "log_lik",
        inc_warmup = FALSE
      )$log_lik
    )

  if (ncol(posterior$beta) > nrow(posterior$beta)) {
    posterior$beta <- t(posterior$beta)
  }
  colnames(posterior$beta) <- colnames(X)
  colnames(posterior$gamma) <- paste0("gamma[", 1:degree, "]")
  colnames(posterior$log_lik) <- paste0("log_lik[", 1:data.n, "]")

  output <- list(
    coefficients = colMeans(posterior$beta),
    bp.param = colMeans(posterior$gamma),
    loglik = colMeans(posterior$log_lik),
    features = features,
    n = data.n,
    nevent = sum(status),
    terms = Terms,
    y = Y,
    formula = formula,
    xlevels = xlevels,
    contrasts = contrasts,
    call = Call,
    means = means,
    sdv = sdv,
    standata = standata,
    posterior = posterior,
    tau_a = 0,
    tau_b = tau,
    stanfit = stanfit
  )

  if (model_flag == "aft") {
    z <- matrix(nrow = nrow(posterior$beta), ncol = length(Y[, 1]))
    for (i in 1:nrow(posterior$beta)) {
      z[i, ] <- (features %*% posterior$beta[i, ])
    }
    output$tau_a <- min(log(Y[, 1])) - max(colMeans(z))
    output$tau_b <- max(log(Y[, 1])) - min(colMeans(z))
  }

  if (null) {
    output$coefficients <- NULL
  }

  output$call <- Call
  output$call$approach <- approach_flag
  output$call$model <- model_flag

  class(output) <- c("spbp")

  return(output)
}
