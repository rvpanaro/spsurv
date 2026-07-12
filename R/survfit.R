#' BP survival, cumulative hazard, SE, and CI at a time grid (internal).
#'
#' @param x Fitted \code{"spbp"} object.
#' @param time Non-negative numeric vector of evaluation times (sorted).
#' @param X Linear predictor design matrix, \code{nrow(X)} prediction curves.
#' @param type CI transform; passed to \code{\link{spsurv:::.survfit_confint}}.
#' @param interval CI level.
#' @param logse Passed to \code{.survfit_confint}.
#' @param interval.type For Bayesian fits only: \code{"hpd"} or \code{"quantile"}.
#' @param monotone Logical; for Bayesian fits, enforce non-increasing band limits over
#'   time (default \code{TRUE} for \code{approach = "bayes"}).
#'
#' @return \code{list(surv, cumhaz, std.err, lower, upper)} matrices \code{length(time)} x \code{nrow(X)}.
#' @importFrom coda mcmc
#' @importFrom stats dbeta pbeta
#' @noRd
.spbp_ph_cumhaz_grad <- function(tau_b, t, gamma, beta, xrow) {
  m <- length(gamma)
  G_row <- vapply(seq_len(m), function(k) pbeta(t / tau_b, k, m - k + 1), numeric(1))
  exp_eta <- if (length(xrow)) exp(sum(xrow * beta)) else 1
  H <- sum(G_row * gamma) * exp_eta
  grad_beta <- if (length(xrow)) H * xrow else numeric(0)
  list(H = H, grad_beta = grad_beta, grad_gamma = G_row * exp_eta)
}

#' @noRd
.spbp_po_R_grad <- function(tau_b, t, gamma, beta, xrow) {
  m <- length(gamma)
  G_row <- vapply(seq_len(m), function(k) pbeta(t / tau_b, k, m - k + 1), numeric(1))
  exp_eta <- if (length(xrow)) exp(sum(xrow * beta)) else 1
  R <- sum(G_row * gamma) * exp_eta
  grad_beta <- if (length(xrow)) R * xrow else numeric(0)
  list(R = R, grad_beta = grad_beta, grad_gamma = G_row * exp_eta)
}

#' @noRd
.spbp_aft_cumhaz_grad <- function(tau_a, tau_b, t, gamma, beta, xrow, P,
                                  X_train = NULL, log_time_train = NULL) {
  m <- length(gamma)
  x_min <- NULL
  x_max <- NULL
  if (!is.null(X_train) && !is.null(log_time_train)) {
    tau_rng <- .spbp_aft_tau_range(beta, X_train, log_time_train, tau_a, tau_b)
    tau_a <- tau_rng$tau_a
    tau_b <- tau_rng$tau_b
    w_train <- log_time_train - as.numeric(X_train %*% beta)
    i_star <- which.min(w_train)[1L]
    j_star <- which.max(w_train)[1L]
    x_min <- X_train[i_star, , drop = TRUE]
    x_max <- X_train[j_star, , drop = TRUE]
  }
  time_aft <- if (length(xrow)) log(t) - sum(xrow * beta) else log(t)
  basis <- .spbp_aft_basis(time_aft = time_aft, tau_a = tau_a, tau_b = tau_b, P = P)
  cumhaz <- as.numeric(basis$G %*% gamma)
  haz <- as.numeric(basis$g %*% gamma)
  if (length(xrow)) {
    if (!is.null(x_min)) {
      range_aft <- tau_b - tau_a
      u <- (time_aft - tau_a) / range_aft
      grad_beta <- haz * (-xrow + x_min - u * (x_min - x_max))
    } else {
      grad_beta <- -haz * xrow
    }
  } else {
    grad_beta <- numeric(0)
  }
  list(
    cumhaz = cumhaz,
    haz = haz,
    grad_beta = grad_beta,
    grad_gamma = as.numeric(basis$G),
    tau_a = tau_a,
    tau_b = tau_b
  )
}

#' @noRd
.spbp_eval_survival <- function(x, time, X, type, interval, logse,
                                interval.type = c("hpd", "quantile"),
                                monotone = FALSE) {
  interval.type <- match.arg(interval.type)
  p <- length(x$coefficients)
  m <- length(x$bp.param)

  if (x$call$approach == "mle") {
    if (p == 0) {
      beta <- 0
    } else {
      beta <- x$coefficients
    }
    gamma <- x$bp.param

    exp_eta <- as.vector(exp(X %*% beta))
    gamma_diag <- .spbp_gamma_information_diagnostics(x)
    if (!gamma_diag$stable) {
      .spbp_warn_gamma_unstable_survfit(x, gamma_diag$kappa_gamma)
    }
    # Full delta-method: propagate joint (beta, gamma) uncertainty.
    # This intentionally includes the gamma part even when p > 0.
    var <- stats::vcov(x, bp.param = TRUE, mask_unstable_gamma = FALSE)

    if (x$call$model == "ph") {
      G <- matrix(sapply(seq_len(m), function(k) pbeta(time / x$tau_b, k, m - k + 1)), nrow = length(time))

      cumhaz <- as.vector(G %*% gamma) %o% exp_eta
      surv <- exp(-cumhaz)
      varH <- matrix(NA_real_, nrow = nrow(cumhaz), ncol = ncol(cumhaz))
      for (i in seq_len(nrow(cumhaz))) {
        for (j in seq_len(ncol(cumhaz))) {
          cg <- .spbp_ph_cumhaz_grad(x$tau_b, time[i], gamma, beta, X[j, , drop = TRUE])
          grad_beta <- cg$grad_beta
          dgamma <- cg$grad_gamma

          if (p > 0) {
            grad <- as.matrix(c(grad_beta, dgamma))
          } else {
            grad <- as.matrix(c(dgamma))
          }

          varH[i, j] <- .spbp_sym_quad_form(var, grad)
        }
      }
      std.err <- sqrt(varH)
    } else if (x$call$model == "po") {
      G <- matrix(sapply(seq_len(m), function(k) pbeta(time / x$tau_b, k, m - k + 1)), nrow = length(time))

      odds <- as.vector(G %*% gamma) %o% exp_eta
      cumhaz <- log(1 + odds)
      surv <- exp(-cumhaz)
      varH <- matrix(NA_real_, nrow = nrow(odds), ncol = ncol(odds))
      for (i in seq_len(nrow(odds))) {
        for (j in seq_len(ncol(odds))) {
          grad_beta <- 1 / (1 + odds[i, j]) * (odds[i, j] * X[j, ])
          dgamma <- 1 / (1 + odds[i, j]) * (G[i, ] * exp_eta[j])

          if (p > 0) {
            grad <- as.matrix(c(grad_beta, dgamma))
          } else {
            grad <- as.matrix(c(dgamma))
          }

          varH[i, j] <- .spbp_sym_quad_form(var, grad)
        }
      }
      std.err <- sqrt(varH)
    } else {
      time_aft <- log(time %o% (1 / exp_eta))
      P <- x$standata$P
      if (is.null(P)) {
        P <- pw.basis(degree = m)
      }
      X_train <- model.matrix(x)
      log_time_train <- log(x$y[, 1])
      tau_rng <- .spbp_aft_tau_range(
        beta = beta,
        X_train = X_train,
        log_time_train = log_time_train,
        tau_a = x$tau_a,
        tau_b = x$tau_b
      )
      tau_a <- tau_rng$tau_a
      tau_b <- tau_rng$tau_b

      cumhaz <- matrix(ncol = nrow(X), nrow = length(time))
      haz <- cumhaz
      g <- list()
      G <- list()

      for (j in seq_len(nrow(X))) {
        basis_j <- .spbp_aft_basis(
          time_aft = time_aft[, j],
          tau_a = tau_a,
          tau_b = tau_b,
          P = P
        )
        g[[j]] <- basis_j$g
        haz[, j] <- as.vector(g[[j]] %*% gamma)
        G[[j]] <- basis_j$G
        cumhaz[, j] <- as.vector(G[[j]] %*% gamma)
      }

      surv <- exp(-cumhaz)
      at_zero <- time == 0
      if (any(at_zero)) {
        surv[at_zero, ] <- 1
        cumhaz[at_zero, ] <- 0
      }
      varH <- matrix(NA_real_, nrow = nrow(cumhaz), ncol = ncol(cumhaz))
      for (i in seq_len(nrow(cumhaz))) {
        for (j in seq_len(ncol(cumhaz))) {
          cg <- .spbp_aft_cumhaz_grad(
            tau_a = tau_a,
            tau_b = tau_b,
            t = time[i],
            gamma = gamma,
            beta = beta,
            xrow = X[j, , drop = TRUE],
            P = P,
            X_train = X_train,
            log_time_train = log_time_train
          )
          grad_beta <- cg$grad_beta
          dgamma <- cg$grad_gamma

          if (p > 0) {
            grad <- as.matrix(c(grad_beta, dgamma))
          } else {
            grad <- as.matrix(c(dgamma))
          }

          varH[i, j] <- .spbp_sym_quad_form(var, grad)
        }
      }
      std.err <- sqrt(varH)
      if (any(at_zero)) {
        std.err[at_zero, ] <- 0
      }
    }

    surv <- surv[, seq_len(nrow(X)), drop = FALSE]
    cumhaz <- cumhaz[, seq_len(nrow(X)), drop = FALSE]
    std.err <- std.err[, seq_len(nrow(X)), drop = FALSE]

    ci <- .survfit_confint(
      p = surv,
      se = std.err,
      logse = logse,
      conf.type = type,
      conf.int = interval
    )
  } else {
    if (p == 0) {
      beta <- array(0, dim = c(1000, 1))
    } else {
      beta <- x$posterior$beta
    }
    gamma <- x$posterior$gamma
    n.samp <- nrow(gamma)
    exp_eta <- matrix(ncol = nrow(X), nrow = n.samp)
    cumhaz <- array(dim = c(length(time), nrow(X), n.samp))
    G <- sapply(seq_len(m), function(k) pbeta(time / x$tau_b, k, m - k + 1))

    if (x$call$model == "ph") {
      for (i in seq_len(n.samp)) {
        exp_eta[i, ] <- as.vector(exp(X %*% beta[i, ]))
        cumhaz[, , i] <- as.vector(G %*% gamma[i, ]) %o% exp_eta[i, ]
      }
    } else if (x$call$model == "po") {
      odds <- cumhaz
      for (i in seq_len(n.samp)) {
        exp_eta[i, ] <- as.vector(exp(X %*% beta[i, ]))
        odds[, , i] <- as.vector(G %*% gamma[i, ]) %o% exp_eta[i, ]
      }
      cumhaz <- log(1 + odds)
    } else {
      time_aft <- array(dim = c(length(time), nrow(X), n.samp))
      P <- x$standata$P
      if (is.null(P)) {
        P <- pw.basis(degree = m)
      }

      for (i in seq_len(n.samp)) {
        exp_eta[i, ] <- as.vector(exp(X %*% beta[i, ]))
        time_aft[, , i] <- log(time %o% (1 / exp_eta[i, ]))
      }

      G <- array(dim = c(length(time), m, nrow(X), n.samp))
      for (j in seq_len(n.samp)) {
        for (i in seq_len(nrow(X))) {
          basis_ij <- .spbp_aft_basis(
            time_aft = time_aft[, i, j],
            tau_a = x$tau_a,
            tau_b = x$tau_b,
            P = P
          )
          G[, , i, j] <- basis_ij$G
          cumhaz[, i, j] <- as.vector(G[, , i, j] %*% gamma[j, ])
        }
      }
    }

    surv <- apply(exp(-cumhaz), c(1, 2), mean)
    std.err <- apply(exp(-cumhaz), c(1, 2), sd)

    surv <- surv[, seq_len(nrow(X)), drop = FALSE]
    std.err <- std.err[, seq_len(nrow(X)), drop = FALSE]

    if (x$call$model == "aft") {
      at_zero <- time == 0
      if (any(at_zero)) {
        surv[at_zero, ] <- 1
        std.err[at_zero, ] <- 0
      }
    }

    surv_draws <- exp(-cumhaz)
    n_time <- dim(surv_draws)[1L]
    n_curve <- dim(surv_draws)[2L]
    ci <- list(
      lower = matrix(NA_real_, nrow = n_time, ncol = n_curve),
      upper = matrix(NA_real_, nrow = n_time, ncol = n_curve)
    )

    .posterior_bounds <- function(surv_draw) {
      surv_draw <- pmax(as.numeric(surv_draw), .Machine$double.eps)
      if (type == "plain") {
        .spbp_posterior_interval(surv_draw, interval, interval.type)
      } else if (type == "log") {
        b <- .spbp_posterior_interval(log(surv_draw), interval, interval.type)
        exp(b)
      } else if (type == "log-log") {
        logh <- -log(surv_draw)
        logh <- pmax(logh, .Machine$double.eps)
        b <- .spbp_posterior_interval(log(logh), interval, interval.type)
        c(exp(-exp(b[2])), exp(-exp(b[1])))
      }
    }

    for (i in seq_len(n_time)) {
      for (j in seq_len(n_curve)) {
        b <- .posterior_bounds(surv_draws[i, j, ])
        ci$lower[i, j] <- b[1]
        ci$upper[i, j] <- b[2]
      }
    }

    if (isTRUE(monotone)) {
      mm <- .spbp_monotone_surv_bands(ci$lower, ci$upper)
      ci$lower <- mm$lower
      ci$upper <- mm$upper
    }

    cumhaz <- apply(cumhaz, c(1, 2), mean)
    cumhaz <- cumhaz[, seq_len(nrow(X)), drop = FALSE]
  }

  list(
    surv = surv,
    cumhaz = cumhaz,
    std.err = std.err,
    lower = ci$lower[, seq_len(nrow(X)), drop = FALSE],
    upper = ci$upper[, seq_len(nrow(X)), drop = FALSE]
  )
}

#' @export
#' @method survfit spbp
#' @title BP-based model survival curves
#' @description Compute survival curves for a fitted \code{\link[spsurv:spbp]{spbp}} model.
#'
#' @param formula An object of class \code{"spbp"} returned by \code{\link[spsurv:spbp]{spbp}}.
#' @param newdata Optional data frame used to obtain survival curves for specific covariate values.
#' @param times Optional evaluation times: a \code{\link[survival:Surv]{Surv}} object (legacy),
#'   or a non-negative \strong{numeric} vector (unique values are used; suitable for smooth
#'   \code{ggplot2::geom_line} plots via \code{\link{predict.spbp}} or \code{\link{as.data.frame.survfitbp}}).
#' @param se.fit Logical; if \code{TRUE}, compute standard errors.
#' @param interval Confidence level for intervals (e.g. \code{0.95}).
#' @param type Character; confidence interval transformation. One of \code{"log"}, \code{"log-log"}, or \code{"plain"}.
#' @param interval.type For Bayesian fits only: \code{"hpd"} (default) or \code{"quantile"} (equal-tailed).
#' @param monotone Logical; for Bayesian fits, enforce non-increasing credible-band limits over
#'   time so ribbons plot smoothly with \code{ggplot2::geom_line}. Defaults to \code{TRUE} when
#'   \code{approach = "bayes"}.
#' @param tidy Logical; if \code{TRUE}, return a \code{data.frame} for ggplot2.
#'   instead of a \code{"survfit"} object.
#' @param ... Further arguments passed to \code{\link[survival:survfit]{survfit.coxph}} for the
#'   reference Cox object (e.g. \code{conf.type} is ignored; use \code{type} instead).
#' @importFrom coda mcmc
#' @importFrom stats model.matrix
#' @return An object of class \code{"survfit"} (with classes \code{survfitbp}, \code{survfitcox}, \code{survfit}).
#' @seealso \code{\link[spsurv:spbp]{spbp}}, \code{\link[survival:survfit]{survfit}},
#'   \code{\link{predict.spbp}}, \code{\link{survfit.spbp}} (\code{as.data.frame} method).
#' @examples
#' library(spsurv)
#' data(veteran, package = "survival")
#' fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran)
#' survfit(fit)
#'
survfit.spbp <- function(formula, newdata = NULL, times = NULL,
                         se.fit = TRUE, interval = .95,
                         type = c("log", "log-log", "plain"),
                         interval.type = c("hpd", "quantile"),
                         monotone = NULL,
                         tidy = FALSE,
                         ...) {
  x <- formula
  type <- match.arg(type)
  interval.type <- match.arg(interval.type)
  if (is.null(monotone)) {
    monotone <- identical(x$call$approach, "bayes")
  }

  data <- .spbp_training_data(x)

  coxfit <- coxph(x$formula, data = data, model = TRUE)
  p <- length(x$coefficients)

  if (is.null(times)) {
    times <- .spbp_default_survfit_times(x)
  }

  newdata_arg <- newdata
  if (is.null(newdata)) {
    km <- survival::survfit(formula = coxfit, ...)
    X <- t(x$means)
  } else {
    km <- survival::survfit(formula = coxfit, newdata = newdata, ...)
    X <- matrix(model.matrix(object = x$formula[-2], xlev = x$xlevels, data = newdata)[, -1], ncol = p)
  }

  if (!is.null(times)) {
    if (inherits(times, "Surv")) {
      km$time <- sort(times[, 1], decreasing = FALSE)
      km$n.event <- times[, 2][order(times[, 1], decreasing = FALSE)]
      km$n.censor <- 1 - km$n.event
      km$n.risk <- length(times[, 1]) - cumsum(km$n.event) + 1
    } else if (is.numeric(times)) {
      tv <- sort(unique(as.numeric(times)))
      tv <- tv[is.finite(tv) & tv >= 0]
      if (!length(tv)) {
        stop("'times' must contain at least one finite non-negative value", call. = FALSE)
      }
      km$time <- tv
      n_t <- length(km$time)
      km$n.event <- rep(0L, n_t)
      km$n.censor <- rep(0L, n_t)
      km$n.risk <- rep(NROW(data), n_t)
    } else {
      stop("'times' must be numeric or a Surv object", call. = FALSE)
    }
  }

  est <- .spbp_eval_survival(
    x,
    time = km$time,
    X = X,
    type = type,
    interval = interval,
    logse = km$logse,
    interval.type = interval.type,
    monotone = monotone
  )

  km$surv <- est$surv
  km$cumhaz <- est$cumhaz
  km$std.err <- est$std.err
  km$lower <- est$lower
  km$upper <- est$upper
  km$conf.type <- type
  km$conf.int <- interval
  km$call <- match.call() ## Call

  if (nrow(X) > 1) {
    colnames(km$surv) <- seq_len(nrow(X))
    colnames(km$cumhaz) <- seq_len(nrow(X))
    colnames(km$std.err) <- seq_len(nrow(X))
    colnames(km$lower) <- seq_len(nrow(X))
    colnames(km$upper) <- seq_len(nrow(X))
  }

  class(km) <- c("survfitbp", "survfitcox", "survfit")

  if (isTRUE(tidy)) {
    df_out <- as.data.frame(km)
    if (!is.null(newdata_arg)) {
      df_out <- .spbp_tidy_survfit_newdata(df_out, newdata_arg)
    }
    return(df_out)
  }

  return(km)
}

#' Merge newdata covariates into tidy survfit output
#' @keywords internal
#' @noRd
.spbp_tidy_survfit_newdata <- function(df, newdata) {
  ids <- unique(df$id)
  n_id <- length(ids)
  n_time <- nrow(df) / n_id
  if (n_time != as.integer(n_time)) {
    return(df)
  }
  idx <- match(as.character(df$id), as.character(seq_len(nrow(newdata))))
  if (any(is.na(idx))) {
    idx <- rep(seq_len(nrow(newdata)), each = n_time)[seq_len(nrow(df))]
  }
  cbind(newdata[idx, , drop = FALSE], df)
}

#' @export
#' @method as.data.frame survfitbp
#' @param x Object from \code{\link{survfit.spbp}}.
#' @param row.names,optional Unused; included for S3 consistency.
#' @return A \code{data.frame} with columns \code{id} (curve index), \code{time}, \code{surv},
#'   \code{lower}, \code{upper}, \code{cumhaz}, \code{std.err}.
#' @describeIn survfit.spbp Tidy survival curves for \code{ggplot2::geom_line} / \code{geom_ribbon}.
#' @examples
#' data(veteran, package = "survival")
#' fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
#' sf <- survfit(fit, times = seq(0, max(veteran$time), length.out = 80))
#' head(as.data.frame(sf))
#'
as.data.frame.survfitbp <- function(x, row.names = NULL, optional = FALSE, ...) {
  tt <- x$time
  surv <- x$surv
  if (!is.matrix(surv)) {
    surv <- matrix(surv, ncol = 1L)
  }
  nc <- ncol(surv)
  ids <- colnames(surv)
  if (is.null(ids)) {
    ids <- as.character(seq_len(nc))
  }

  lower <- x$lower
  upper <- x$upper
  if (!is.matrix(lower)) {
    lower <- matrix(lower, ncol = nc)
  }
  if (!is.matrix(upper)) {
    upper <- matrix(upper, ncol = nc)
  }

  ch <- x$cumhaz
  se <- x$std.err
  if (!is.matrix(ch)) {
    ch <- matrix(ch, ncol = nc)
  }
  if (!is.matrix(se)) {
    se <- matrix(se, ncol = nc)
  }

  data.frame(
    id = rep(ids, each = length(tt)),
    time = rep(tt, times = nc),
    surv = as.vector(surv),
    lower = as.vector(lower),
    upper = as.vector(upper),
    cumhaz = as.vector(ch),
    std.err = as.vector(se),
    stringsAsFactors = FALSE
  )
}

#' @export
#' @method predict spbp
#' @title Predicted survival as a tidy data frame
#' @description Survival (and optional CI) on a time grid, as a \code{data.frame} for
#'   \code{ggplot2::geom_line}.
#'
#' @param object Fitted \code{"spbp"} from \code{\link{spbp}} / \code{\link{bpph}} / etc.
#' @param newdata Optional \code{data.frame} of covariate profiles (same convention as \code{\link{survfit.spbp}}).
#' @param times Time grid. Default: a dense sequence from 0 to the maximum observed
#'   time (suitable for smooth \code{ggplot2::geom_line} / \code{geom_ribbon} plots).
#'   Pass a long \code{seq(...)} to override resolution; use observed event times only
#'   if stepwise curves are intended.
#' @param interval,type Passed to \code{\link{survfit.spbp}}.
#' @param interval.type,monotone Passed to \code{\link{survfit.spbp}} (Bayesian fits).
#' @param ... Passed to \code{\link{survfit.spbp}}.
#' @return Same structure as \code{\link{as.data.frame.survfitbp}}.
#' @seealso \code{\link{survfit.spbp}}
#' @examples
#' data(veteran, package = "survival")
#' fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
#' pr <- predict(fit, times = seq(0, 400, by = 2))
#' \dontrun{
#'   ggplot2::ggplot(pr, ggplot2::aes(time, surv)) + ggplot2::geom_line()
#' }
#'
predict.spbp <- function(object, newdata = NULL, times = NULL, interval = .95,
                         type = c("log", "log-log", "plain"),
                         interval.type = c("hpd", "quantile"),
                         monotone = NULL,
                         ...) {
  type <- match.arg(type)
  interval.type <- match.arg(interval.type)
  if (is.null(times)) {
    times <- .spbp_default_survfit_times(object)
  } else {
    times <- sort(unique(as.numeric(times)))
    times <- times[is.finite(times) & times >= 0]
    if (!length(times)) {
      stop("'times' must contain at least one finite non-negative value", call. = FALSE)
    }
  }

  if (is.null(newdata)) {
    survfit(
      object,
      times = times,
      interval = interval,
      type = type,
      interval.type = interval.type,
      monotone = monotone,
      tidy = TRUE,
      ...
    )
  } else {
    survfit(
      object,
      newdata = newdata,
      times = times,
      interval = interval,
      type = type,
      interval.type = interval.type,
      monotone = monotone,
      tidy = TRUE,
      ...
    )
  }
}

#' @export
#' @method residuals spbp
#' @title BP based models residuals.
#' @description Residuals for a fitted \code{\link[spsurv]{spbp}} model.
#' @param object an object of class `spbp` result of a \code{\link[spsurv]{spbp}} fit.
#' @param type type of residuals, default is "cox-snell"
#' @seealso \code{\link[spsurv:spbp]{spbp}}, \code{\link[survival:survfit]{spbp}}.
#' @param ... arguments passed to parent method.
#' @examples
#'
#' library("spsurv")
#' data("veteran", package = "survival")
#'
#' fit <- bpph(Surv(time, status) ~ karno + factor(celltype),
#'   data = veteran
#' )
#'
#' residuals(fit)
residuals.spbp <- function(object, type = c("martingale", "deviance", "cox-snell", "coxsnell"), ...) {
  type <- match.arg(type)
  if (identical(type, "coxsnell")) {
    type <- "cox-snell"
  }

  p <- length(object$coefficients)

  if (p > 0) {
    X <- model.matrix(object)
  } else {
    X <- t(object$means)
  }

  m <- length(object$bp.param)
  cumhaz <- NULL

  if (object$call$approach == "mle") {
    if (p == 0) {
      beta <- 0
    } else {
      beta <- object$coefficients
    }

    gamma <- object$bp.param
    exp_eta <- as.vector(exp(X %*% beta))

    if (object$call$model == "ph") {
      G <- matrix(sapply(seq_len(m), function(k) pbeta(object$y[, 1] / object$tau_b, k, m - k + 1)), nrow = length(object$y[, 1]))
      cumhaz <- as.vector(G %*% gamma) * exp_eta
    } else if (object$call$model == "po") {
      G <- matrix(sapply(seq_len(m), function(k) pbeta(object$y[, 1] / object$tau_b, k, m - k + 1)), nrow = length(object$y[, 1]))
      odds <- as.vector(G %*% gamma) * exp_eta
      cumhaz <- log(1 + odds)
    } else {
      n_obs <- length(object$y[, 1])
      P <- object$standata$P
      if (is.null(P)) {
        P <- pw.basis(degree = m)
      }
      lp <- if (p > 0) {
        as.vector(X %*% as.matrix(beta, ncol = 1L))
      } else {
        rep(0, n_obs)
      }
      time <- log(object$y[, 1]) - lp
      G <- .spbp_aft_basis(
        time_aft = time,
        tau_a = object$tau_a,
        tau_b = object$tau_b,
        P = P
      )$G
      cumhaz <- as.vector(G %*% gamma)
    }
  } else {
    beta <- object$posterior$beta
    gamma <- object$posterior$gamma
    n.samp <- nrow(beta)
    exp_eta <- matrix(nrow = n.samp, ncol = length(object$y[, 1]))
    cumhaz <- array(dim = c(n.samp, length(object$y[, 1])))
    G <- sapply(seq_len(m), function(k) pbeta(object$y[, 1] / object$tau_b, k, m - k + 1))

    if (object$call$model == "ph") {
      for (i in seq_len(n.samp)) {
        exp_eta[i, ] <- as.vector(exp(X %*% beta[i, ]))
        cumhaz[i, ] <- as.vector(G %*% gamma[i, ]) * exp_eta[i, ]
      }
    } else if (object$call$model == "po") {
      odds <- cumhaz
      for (i in seq_len(n.samp)) {
        exp_eta[i, ] <- as.vector(exp(X %*% beta[i, ]))
        odds[i, ] <- as.vector(G %*% gamma[i, ]) * exp_eta[i, ]
      }
      cumhaz <- log(1 + odds)
    } else {
      time <- array(dim = c(n.samp, length(object$y[, 1])))
      P <- object$standata$P
      if (is.null(P)) {
        P <- pw.basis(degree = m)
      }

      for (i in seq_len(n.samp)) {
        exp_eta[i, ] <- as.vector(exp(X %*% beta[i, ]))
        time[i, ] <- log(object$y[, 1] / exp_eta[i, ])
      }

      G <- array(dim = c(n.samp, length(object$y[, 1]), m))
      for (i in seq_len(n.samp)) {
        G[i, , ] <- .spbp_aft_basis(
          time_aft = time[i, ],
          tau_a = object$tau_a,
          tau_b = object$tau_b,
          P = P
        )$G
        cumhaz[i, ] <- as.vector(G[i, , ] %*% gamma[i, ])
      }
    }
    cumhaz <- apply(cumhaz, 2, mean)
  }

  cumhaz <- as.vector(cumhaz)
  names(cumhaz) <- names(object$y[, 1])
  delta <- object$y[, 2]

  if (type == "cox-snell") {
    return(cumhaz)
  } else if (type == "martingale") {
    return(delta - cumhaz)
  } else {
    m <- delta - cumhaz
    sign(m) * sqrt(-2 * (m + delta * log(delta - m)))
  }
}
