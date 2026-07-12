#' Semiparametric Survival Analysis Using Bernstein Polynomial
#'
#' Fits Bernstein PH, PO, or AFT models to survival data via Stan (MLE or
#' Bayesian).
#'
#' @title spbp: The BP Based Survival Analysis Function
#' @param formula a \code{\link[survival]{Surv}} response with event time,
#'   censoring status, and optional covariates.
#' @param ... Arguments passed to \code{\link{spbp.default}}, including
#'   \code{data}, \code{model}, \code{approach}, and \code{degree}. Further
#'   arguments in \code{...} are passed to \code{rstan::optimizing} (MLE) or
#'   \code{rstan::sampling} (Bayes), e.g. \code{iter}, \code{chains}, \code{init}.
#' @details
#' The generic dispatches to \code{\link{spbp.default}} for formula objects.
#' Convenience wrappers \code{\link{bpph}}, \code{\link{bppo}}, and
#' \code{\link{bpaft}} fix the model family. See
#' \code{vignette("getting-started", package = "spsurv")} for a tutorial,
#' \code{vignette("model-families", package = "spsurv")} for PH / PO / AFT
#' comparison, and \code{vignette("bp-degree", package = "spsurv")} for
#' choosing the Bernstein polynomial degree.
#' @return An object of class \code{"spbp"}. See \code{\link{spbp.default}} for
#'   the list of components (\code{coefficients}, \code{bp.param},
#'   \code{degree}, etc.).
#' @seealso \code{\link{spbp.default}}, \code{\link{bpph}}, \code{\link{bppo}},
#'   \code{\link{bpaft}}, \code{\link{bernstein}}
#' @examples
#'
#' library("spsurv")
#' data("veteran", package = "survival")
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
#' @importFrom stats terms

spbp <- function(formula, ...) {
  UseMethod("spbp", formula)
}

#' @title spbp: The BP Based Semiparametric Survival Analysis Function
#' @param formula a Surv object with time to event, status and explanatory terms
#' @param degree Bernstein polynomial degree (integer). If omitted and neither
#'   \code{dist} nor \code{baseline} supplies \code{\link{bernstein}(m)}, the
#'   default is \code{ceiling(sqrt(n))} where \code{n} is the number of rows in
#'   \code{data}.
#' @param data a data.frame object
#' @param approach Bayesian or Maximum Likelihood estimation methods; default is \code{"mle"}
#' @param model Bernstein PH (\code{"ph"}), PO (\code{"po"}), or AFT (\code{"aft"}) model; default is \code{"ph"}
#' @param priors prior settings for the Bayesian approach; `normal` or `cauchy` for beta; `lognormal` or `loglogistic` for gamma (BP coefficients)
#' @param scale logical; indicates whether to center and scale the data
#' @param dist optional baseline specification; use \code{\link{bernstein}(m)} for the Bernstein polynomial degree
#' @param baseline optional alias for \code{dist}
#' @param ... further arguments passed to \code{rstan::optimizing} (MLE) or
#'   \code{rstan::sampling} (Bayes), e.g. \code{iter}, \code{warmup}, \code{init}.
#' @param cores number of core threads to use (Bayes sampling)
#' @param verbose passed to Stan
#' @param chains number of MCMC chains (Bayes)
#' @details
#' Right-censored survival data are modeled with a Bernstein-polynomial baseline
#' and regression on covariates. With \code{approach = "mle"}, parameters are
#' estimated by Stan's optimizer and approximate inference uses the Hessian when
#' available. With \code{approach = "bayes"}, posterior samples are drawn with
#' NUTS; use \code{summary}, \code{\link{tidy.spbp}}, and \code{\link{glance.spbp}}
#' for output. Covariates are centered and scaled when \code{scale = TRUE}
#' (default).
#'
#' The returned object includes:
#' \describe{
#'   \item{\code{coefficients}}{Regression estimates on the original covariate scale.}
#'   \item{\code{bp.param}}{Bernstein baseline coefficients (\code{gamma}).}
#'   \item{\code{degree}}{Polynomial degree used (also in \code{call$degree}).}
#'   \item{\code{loglik}}{MLE: intercept-only and full-model log-likelihoods;
#'     Bayes: posterior mean pointwise log-likelihoods.}
#'   \item{\code{call}}{Matched call with \code{approach}, \code{model}, and \code{degree}.}
#' }
#' @return An object of class \code{spbp}. Component \code{degree} records the
#'   Bernstein polynomial degree used in the fit (also stored in \code{call$degree}).
#' @seealso \code{\link{bpph}}, \code{\link{bppo}}, \code{\link{bpaft}},
#'   \code{\link{bernstein}}, \code{\link{summary.spbp}}
#' @method spbp default
#' @export
#' @importFrom rstan stan sampling optimizing
#' @importFrom survival Surv frailty
#' @importFrom coda HPDinterval
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
           cores = .spbp_default_cores(),
           scale = TRUE,
           dist = NULL,
           baseline = NULL,
           verbose = FALSE, chains = 4,
           ...) {
    # ---------------Definitions + error handling  ---------------
    Call <- match.call()

    if (anyNA(cores) || length(cores) < 1L || cores[1L] < 1L) {
      dc <- suppressWarnings(parallel::detectCores())
      cores <- if (length(dc) != 1L || is.na(dc) || dc < 2L) {
        1L
      } else {
        min(as.integer(dc - 1L), 4L)
      }
    }
    cores <- as.integer(cores[1L])

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

    aux <- match(c("formula", "data"), names(Call), nomatch = 0)
    temp <- Call[c(1, aux)] # keep important args
    temp[[1L]] <- quote(stats::model.frame) # model frame call
    special <- c("frailty", "frailty.gamma", "frailty.gaussian", "frailty.t")
    temp$formula <- terms(formula, special, data = data)
    stanArgs <- list(...)

    mf <- eval(temp, parent.frame())
    Terms <- terms(mf)
    Y <- model.extract(mf, "response")
    if (!inherits(Y, "Surv")) {
      stop("Response must be a survival object")
    }
    type <- attr(Y, "type")

    time <- as.vector(Y[, 1])
    tau <- max(time)
    status <- as.vector(Y[, 2])

    degree <- .spbp_resolve_degree(degree, dist = dist, baseline = baseline)
    if (is.null(degree)) {
      degree <- ceiling(nrow(data)^(0.5))
    }

    # Frailty: id, distribution, column index
    id <- NULL
    if (!is.null(attr(temp$formula, "specials")$frailty)) {
      frailty_idx <- attr(temp$formula, "specials")$frailty
      id <- model.matrix(formula)[, attr(temp$formula, "specials")$frailty]
      rand <- 1L
    } else if (!is.null(attr(temp$formula, "specials")$frailty.gamma)) {
      frailty_idx <- attr(temp$formula, "specials")$frailty.gamma
      id <- model.matrix(formula)[, attr(temp$formula, "specials")$frailty.gamma]
      rand <- 1L
    } else if (!is.null(attr(temp$formula, "specials")$frailty.gaussian)) {
      frailty_idx <- attr(temp$formula, "specials")$frailty.gaussian
      id <- model.matrix(formula)[, attr(temp$formula, "specials")$frailty.gaussian]
      rand <- 2L
    } else if (!is.null(attr(temp$formula, "specials")$frailty.t)) {
      frailty_idx <- attr(temp$formula, "specials")$frailty.t
      id <- model.matrix(formula)[, attr(temp$formula, "specials")$frailty.t]
      rand <- 3L
    } else {
      rand <- 0L
      frailty_idx <- NULL
    }

    # Parse priors (beta, gamma, frailty)
    if (length(priors$beta) > 0) {
      betap <- lapply(priors$beta, read_prior)
    } else {
      betap <- list(c("normal", "0", "5"))
    }
    priordist_beta <- sapply(betap, `[[`, 1)
    location_beta <- sapply(betap, `[[`, 2)
    scale_beta <- sapply(betap, `[[`, 3)
    if (length(priors$gamma) > 0) {
      gammap <- lapply(priors$gamma, read_prior)
    } else {
      gammap <- list(c("lognormal", "0", "5"))
    }
    priordist_gamma <- sapply(gammap, `[[`, 1)
    location_gamma <- sapply(gammap, `[[`, 2)
    scale_gamma <- sapply(gammap, `[[`, 3)
    if (length(priors$frailty) > 0) {
      frailtyp <- lapply(priors$frailty, read_prior)
    } else {
      frailtyp <- list(c("gamma", "1", "1"))
    }
    priordist_frailty <- sapply(frailtyp, `[[`, 1)
    par1_frailty <- sapply(frailtyp, `[[`, 2)
    par2_frailty <- sapply(frailtyp, `[[`, 3)

    # Model frame / response validation
    if (nrow(mf) == 0) stop("Only missing observations")
    if (!type %in% c("right", "counting")) {
      stop(paste0("spsurv doesn't support \"", type, "\" survival data"))
    }
    if (length(attr(Terms, "variables")) > 2) {
      ytemp <- terms.inner(formula)[1:2]
      xtemp <- terms.inner(formula)[-c(1, 2)]
      if (any(!is.na(match(xtemp, ytemp)))) {
        stop("a variable appears on both the left and right sides of the formula")
      }
    }

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

    if (scale && null == 0L) {
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
        expr = .spbp_mle(standata = standata, hessian = TRUE, verbose = verbose, ...),
        error = function(e) {
          warning(e)
          return(NaN)
        }
      )
    } else {
      tryCatch(
        expr = .spbp_bayes(standata = standata, hessian = FALSE, verbose = verbose, chains = chains, ...),
        error = function(e) {
          warning(e)
          return(NaN)
        }
      )
  }
}

#' Initial values from prior for MLE optimizer (internal)
#' @param standata Stan data list with p, m, location/scale and priordist vectors.
#' @return List with elements \code{beta} and \code{gamma}.
#' @keywords internal
#' @noRd
.spbp_initial_values <- function(standata) {
  p <- standata$p
  m <- standata$m

  beta_init <- numeric(p)
  for (j in seq_len(p)) {
    # normal(location, scale): one draw from prior (in 95% interval w.h.p.)
    beta_init[j] <- stats::rnorm(1, mean = 0, sd = 1)
  }

  gamma_init <- numeric(m)
  for (j in seq_len(m)) {
    gamma_init[j] <- exp(stats::rnorm(1, mean = 0, sd = 1))
  }

  # rstan::optimizing expects a *flat* list (beta[1], ..., gamma[1], ...), not
  # list(beta = ..., gamma = ...); the latter breaks for vector length 1
  # (declared dim 1, found empty) in recent rstan.
  init_beta <- if (p > 0L) {
    stats::setNames(as.list(beta_init), paste0("beta[", seq_len(p), "]"))
  } else {
    list()
  }
  init_gamma <- stats::setNames(as.list(gamma_init), paste0("gamma[", seq_len(m), "]"))
  c(init_beta, init_gamma)
}

#' MLE fit for spbp (internal)
#' @keywords internal
#' @noRd
.spbp_mle <-
  function(standata, hessian = TRUE, verbose = FALSE, ...) {
    e <- parent.frame() # "sourcing" the parent.frame
    vnames <- objects(, envir = e) # variable names in parent frame

    for (n in vnames) assign(n, get(n, e))

    if (!is.null(frailty_idx)) {
      standata$X <- X[, -frailty_idx]
      message("Frailty ignored, change approach to `bayes` for frailty estimation.")
    }

    dots <- list(...)
    user_init <- dots$init
    dots$init <- NULL

    hessian_run <- hessian
    attempt <- 0L
    stanfit <- list(return_code = 70)

    fit_opt <- function(init, h) {
      suppressWarnings(do.call(rstan::optimizing, c(list(
        object  = stanmodels$spbp,
        data    = standata,
        init    = init,
        hessian = h,
        verbose = verbose
      ), dots)))
    }

    while (stanfit$return_code != 0 && attempt < 3L) {
      attempt <- attempt + 1L
      # Initial values: caller may pass init=... to rstan::optimizing (e.g. init = 0 in tests);
      # otherwise draw a flat list (beta[1], ..., gamma[1], ...) required by rstan::optimizing.
      init <- if (!is.null(user_init)) {
        user_init
      } else {
        .spbp_initial_values(standata)
      }
      stanfit <- tryCatch(
        fit_opt(init, hessian_run),
        error = function(e) {
          msg <- conditionMessage(e)
          if (isTRUE(hessian_run) && isTRUE(hessian) &&
                grepl("optimHess|non-finite|hessian|finite difference", msg, ignore.case = TRUE)) {
            hessian_run <<- FALSE
            warning(
              "Hessian computation failed; refitting without Hessian. ",
              call. = FALSE
            )
            fit_opt(init, FALSE)
          } else {
            stop(e)
          }
        }
      )
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

    n_full <- 1L + standata$p + standata$m
    n_bg <- standata$p + standata$m
    h_mat <- stanfit$hessian
    has_hess <- !is.null(h_mat) && is.matrix(h_mat) && nrow(h_mat) == ncol(h_mat)
    dim_ok <- has_hess && (nrow(h_mat) == n_full || nrow(h_mat) == n_bg)
    if (!hessian_run || !dim_ok) {
      stanfit$hessian <- matrix(NA_real_, n_full, n_full)
    }

    nulldata <- standata
    nulldata$X <- matrix(0, ncol = ncol(X), nrow = data.n)

    nullfit <- tryCatch(
      do.call(
        rstan::optimizing,
        c(
          list(
            object = stanmodels$spbp,
            data = nulldata,
            hessian = hessian_run,
            init = if (!is.null(user_init)) user_init else "random"
          ),
          dots
        )
      ),
      error = function(e) {
        return(list(value = NULL))
      }
    )

    output <- list(
      coefficients = coef,
      bp.param = gamma,
      degree = degree,
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
      stanfit = stanfit,
      data = data
    )

    if (model_flag == "aft") {
      output$tau_a <- min(log(Y[, 1]) - features %*% coef)
      output$tau_b <- max(log(Y[, 1]) - features %*% coef)
    }

    if (null) {
      output$coefficients <- NULL
      output$hessian <- stanfit$hessian[-1, -1]
      output$features <- NULL
    }

    output$call$approach <- approach_flag
    output$call$model <- model_flag
    output$call$degree <- degree

    class(output) <- c("spbp")

    if (isTRUE(verbose)) {
      message("Priors are ignored because the MLE approach is used.")
    }

    if (isTRUE(getOption("spsurv.store_residuals", TRUE))) {
      output$residuals <- residuals(output)
    }

    return(output)
  }

#' Bayesian fit for spbp (internal)
#'
#' @keywords internal
.spbp_bayes <- function(standata, hessian = FALSE, verbose = FALSE, chains = 1, ...) {
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
  posterior$draw_indices <- .spbp_posterior_draw_indices(stanfit, nrow(posterior$beta))

  output <- list(
    coefficients = colMeans(posterior$beta),
    bp.param = colMeans(posterior$gamma),
    degree = degree,
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
    stanfit = stanfit,
    data = data
  )

  if (model_flag == "aft") {
    z <- matrix(nrow = nrow(posterior$beta), ncol = length(Y[, 1]))
    for (i in seq_len(nrow(posterior$beta))) {
      z[i, ] <- (features %*% posterior$beta[i, ])
    }
    output$tau_a <- min(log(Y[, 1]) - colMeans(z))
    output$tau_b <- max(log(Y[, 1]) - colMeans(z))
  }

  if (null) {
    output$coefficients <- NULL
  }

  output$call <- Call
  output$call$approach <- approach_flag
  output$call$model <- model_flag
  output$call$degree <- degree

  class(output) <- c("spbp")

  return(output)
}
