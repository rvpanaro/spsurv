#' Convert spbp Bayesian fit to posterior draws
#'
#' @param x A fitted \code{"spbp"} object from \code{approach = "bayes"}.
#' @param variables Character vector of posterior components to include
#'   (\code{"beta"}, \code{"gamma"}, \code{"log_lik"}).
#' @param ... Not used.
#' @return A \code{posterior::draws_df} object on the back-transformed scale.
#' @export
as_draws_df.spbp <- function(x, variables = c("beta", "gamma"), ...) {
  if (!identical(x$call$approach, "bayes") || is.null(x$posterior)) {
    stop("as_draws_df() requires a Bayesian spbp fit.", call. = FALSE)
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required.", call. = FALSE)
  }
  variables <- match.arg(variables, c("beta", "gamma", "log_lik"), several.ok = TRUE)
  idx <- x$posterior$draw_indices
  if (is.null(idx)) {
    idx <- data.frame(
      .chain = 1L,
      .iteration = seq_len(nrow(x$posterior$beta)),
      .draw = seq_len(nrow(x$posterior$beta)),
      stringsAsFactors = FALSE
    )
  }
  out <- idx
  if ("beta" %in% variables && !is.null(x$posterior$beta)) {
    beta_df <- as.data.frame(x$posterior$beta)
    names(beta_df) <- paste0("beta[", names(beta_df), "]")
    out <- cbind(out, beta_df)
  }
  if ("gamma" %in% variables && !is.null(x$posterior$gamma)) {
    gamma_df <- as.data.frame(x$posterior$gamma)
    names(gamma_df) <- colnames(x$posterior$gamma)
    out <- cbind(out, gamma_df)
  }
  if ("log_lik" %in% variables && !is.null(x$posterior$log_lik)) {
    ll_df <- as.data.frame(x$posterior$log_lik)
    names(ll_df) <- colnames(x$posterior$log_lik)
    out <- cbind(out, ll_df)
  }
  posterior::as_draws_df(out)
}

#' Tidy posterior draws for spbp fits
#'
#' @param object A Bayesian \code{"spbp"} fit.
#' @param ... Passed to \code{\link{as_draws_df.spbp}}.
#' @return A \code{posterior::draws_df} object.
#' @keywords internal
#' @noRd
tidy_draws.spbp <- function(object, ...) {
  as_draws_df.spbp(object, ...)
}

#' Spread spbp draws into long format
#'
#' @param object A Bayesian \code{"spbp"} fit.
#' @param ... Expressions passed to \code{tidybayes::spread_draws}.
#' @return A tibble.
#' @keywords internal
#' @noRd
spread_draws.spbp <- function(object, ...) {
  if (!requireNamespace("tidybayes", quietly = TRUE)) {
    stop("Package 'tidybayes' is required.", call. = FALSE)
  }
  draws <- as_draws_df.spbp(object)
  tidybayes::spread_draws(draws, ...)
}

#' Gather spbp draws into long format
#'
#' @param object A Bayesian \code{"spbp"} fit.
#' @param ... Expressions passed to \code{tidybayes::gather_draws}.
#' @return A tibble.
#' @keywords internal
#' @noRd
gather_draws.spbp <- function(object, ...) {
  if (!requireNamespace("tidybayes", quietly = TRUE)) {
    stop("Package 'tidybayes' is required.", call. = FALSE)
  }
  draws <- as_draws_df.spbp(object)
  tidybayes::gather_draws(draws, ...)
}

#' @keywords internal
#' @noRd
.spbp_bayes_surv_draws_array <- function(x, time, X) {
  if (!identical(x$call$approach, "bayes") || is.null(x$posterior)) {
    stop("Posterior survival draws require a Bayesian spbp fit.", call. = FALSE)
  }
  p <- length(x$coefficients)
  m <- length(x$bp.param)
  beta <- x$posterior$beta
  gamma <- x$posterior$gamma
  n.samp <- nrow(gamma)
  cumhaz <- array(dim = c(length(time), nrow(X), n.samp))
  G <- sapply(seq_len(m), function(k) pbeta(time / x$tau_b, k, m - k + 1))

  if (x$call$model == "ph") {
    for (i in seq_len(n.samp)) {
      exp_eta <- as.vector(exp(X %*% beta[i, ]))
      cumhaz[, , i] <- as.vector(G %*% gamma[i, ]) %o% exp_eta
    }
  } else if (x$call$model == "po") {
    odds <- cumhaz
    for (i in seq_len(n.samp)) {
      exp_eta <- as.vector(exp(X %*% beta[i, ]))
      odds[, , i] <- as.vector(G %*% gamma[i, ]) %o% exp_eta
    }
    cumhaz <- log(1 + odds)
  } else {
    time_aft <- array(dim = c(length(time), nrow(X), n.samp))
    P <- x$standata$P
    if (is.null(P)) {
      P <- pw.basis(degree = m)
    }
    for (i in seq_len(n.samp)) {
      exp_eta <- as.vector(exp(X %*% beta[i, ]))
      time_aft[, , i] <- log(time %o% (1 / exp_eta))
    }
    for (j in seq_len(n.samp)) {
      for (i in seq_len(nrow(X))) {
        basis_ij <- .spbp_aft_basis(
          time_aft = time_aft[, i, j],
          tau_a = x$tau_a,
          tau_b = x$tau_b,
          P = P
        )
        cumhaz[, i, j] <- as.vector(basis_ij$G %*% gamma[j, ])
      }
    }
  }

  exp(-cumhaz)
}

#' Posterior survival curves in long tidy format
#'
#' @param object A Bayesian \code{"spbp"} fit.
#' @param times Numeric vector of evaluation times.
#' @param newdata Optional covariate profiles.
#' @param ... Not used.
#' @return A \code{data.frame} with \code{.chain}, \code{.iteration}, \code{.draw},
#'   \code{time}, \code{surv}, and optional covariate columns.
#' @export
spread_surv_draws.spbp <- function(object, times, newdata = NULL, ...) {
  if (!identical(object$call$approach, "bayes")) {
    stop("spread_surv_draws() requires approach = 'bayes'.", call. = FALSE)
  }
  times <- sort(unique(as.numeric(times)))
  times <- times[is.finite(times) & times >= 0]
  if (!length(times)) {
    stop("'times' must contain at least one finite non-negative value.", call. = FALSE)
  }
  if (is.null(newdata)) {
    p <- length(object$coefficients)
    if (p == 0L) {
      X <- matrix(numeric(0), nrow = 1L, ncol = 0L)
    } else {
      X <- t(object$means)
    }
    newdata <- NULL
  } else {
    X <- .spbp_newdata_matrix(object, newdata)
  }

  surv_draws <- .spbp_bayes_surv_draws_array(object, time = times, X = X)
  idx <- object$posterior$draw_indices
  if (is.null(idx)) {
    n_draws <- dim(surv_draws)[3L]
    idx <- data.frame(
      .chain = 1L,
      .iteration = seq_len(n_draws),
      .draw = seq_len(n_draws),
      stringsAsFactors = FALSE
    )
  }

  n_time <- length(times)
  n_id <- nrow(X)
  n_draws <- dim(surv_draws)[3L]
  rows <- n_time * n_id * n_draws
  out <- data.frame(
    .chain = rep(idx$.chain, each = n_time * n_id),
    .iteration = rep(idx$.iteration, each = n_time * n_id),
    .draw = rep(idx$.draw, each = n_time * n_id),
    time = rep(rep(times, each = n_id), times = n_draws),
    surv = as.vector(surv_draws),
    id = rep(rep(seq_len(n_id), n_time), times = n_draws),
    stringsAsFactors = FALSE
  )
  if (!is.null(newdata)) {
    cov_cols <- newdata[rep(seq_len(n_id), n_time * n_draws), , drop = FALSE]
    out <- cbind(cov_cols, out)
  }
  out
}
