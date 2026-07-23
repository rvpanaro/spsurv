#'
#' @export
terms.inner <- function(x, ...) {
  if (is(x, "formula")) {
    c(terms.inner(x[[2]]), terms.inner(x[[3]]))
  } else if (is(x, "call") && (x[[1]] != as.name("$") &&
    x[[1]] != as.name("["))) {
    if (x[[1]] == "+" || x[[1]] == "*" || x[[1]] == "-") {
      c(terms.inner(x[[2]]), terms.inner(x[[3]]))
    } else if (x[[1]] == as.name("Surv") || x[[1]] == as.name("rand")) {
      unlist(lapply(x[-1], terms.inner))
    } else {
      terms.inner(x[[2]])
    }
  } else {
    (deparse(x))
  }
}

read_prior <- function(prior) {
  aux <- unlist(strsplit(prior, "\\("))
  dist <- aux[1]
  aux2 <- unlist(strsplit(aux[2], "\\)"))[1]
  val <- unlist(strsplit(aux2, ","))
  c(dist, val)
}

#' Training data for \code{survfit.spbp} and similar (Cox refit).
#'
#' Prefers \code{x$data} (the data frame passed to \code{spbp}); otherwise
#' evaluates \code{call$data} in \code{environment(terms)}, then walks the call
#' stack so symbols like \code{data = veteran2} resolve after \code{readRDS}
#' when the terms environment is empty.
#'
#' @keywords internal
#' @noRd
.spbp_training_data <- function(x) {
  if (!is.null(x$data)) {
    return(x$data)
  }
  cd <- x$call$data
  env_terms <- environment(x$terms)
  if (!is.null(env_terms)) {
    val <- tryCatch(eval(cd, envir = env_terms), error = function(e) NULL)
    if (!is.null(val)) {
      return(val)
    }
  }
  for (i in sys.nframe():1) {
    fr <- sys.frame(i)
    val <- tryCatch(eval(cd, envir = fr), error = function(e) NULL)
    if (!is.null(val) && inherits(val, "data.frame")) {
      return(val)
    }
  }
  if (!is.null(env_terms)) {
    return(eval(cd, envir = env_terms))
  }
  stop(
    "unable to recover training data for this 'spbp' object; refit the model or upgrade spsurv so fits store the training data.",
    call. = FALSE
  )
}

#' Default MCMC core count for \code{spbp} (avoids \code{NA} from \code{detectCores}).
#' @keywords internal
#' @noRd
.spbp_default_cores <- function() {
  dc <- suppressWarnings(parallel::detectCores())
  if (length(dc) != 1L || is.na(dc) || dc < 2L) {
    1L
  } else {
    min(as.integer(dc - 1L), 4L)
  }
}

#' Default evaluation grid for smooth survival curves (internal).
#' @keywords internal
#' @noRd
.spbp_survfit_train_nevent <- function(x, km) {
  if (!is.null(x$nevent)) {
    return(as.integer(x$nevent))
  }
  if (!is.null(km$n.event)) {
    return(as.integer(sum(km$n.event, na.rm = TRUE)))
  }
  NA_integer_
}

#' @keywords internal
#' @noRd
.spbp_default_survfit_times <- function(x, length.out = NULL) {
  t_max <- max(x$y[, 1L], na.rm = TRUE)
  if (!is.finite(t_max) || t_max <= 0) {
    return(0)
  }
  if (is.null(length.out)) {
    n <- x$n
    length.out <- max(100L, min(500L, as.integer(round(10 * sqrt(n)))))
  }
  seq(0, t_max, length.out = length.out)
}

#' Enforce non-increasing pointwise survival band limits over time (internal).
#' @keywords internal
#' @noRd
.spbp_monotone_surv_bands <- function(lower, upper) {
  if (!is.matrix(lower)) {
    lower <- matrix(lower, ncol = 1L)
  }
  if (!is.matrix(upper)) {
    upper <- matrix(upper, ncol = 1L)
  }
  for (j in seq_len(ncol(lower))) {
    lower[, j] <- cummin(lower[, j])
    upper[, j] <- cummin(upper[, j])
    upper[, j] <- pmax(upper[, j], lower[, j])
  }
  list(
    lower = pmax(lower, 0),
    upper = pmin(pmax(upper, lower), 1)
  )
}

#' Draw index metadata for Bayesian spbp fits
#' @keywords internal
#' @noRd
.spbp_posterior_draw_indices <- function(stanfit, n_draws) {
  n_chains <- stanfit@sim$chains
  n_iter <- stanfit@sim$iter - stanfit@sim$warmup
  if (n_chains >= 1L && n_iter >= 1L && n_chains * n_iter == n_draws) {
    chain <- rep(seq_len(n_chains), each = n_iter)
    iteration <- rep(seq_len(n_iter), times = n_chains)
    return(data.frame(
      .chain = chain,
      .iteration = iteration,
      .draw = seq_len(n_draws),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    .chain = 1L,
    .iteration = seq_len(n_draws),
    .draw = seq_len(n_draws),
    stringsAsFactors = FALSE
  )
}

#' Posterior interval for one vector of draws (internal).
#' @keywords internal
#' @noRd
.spbp_posterior_interval <- function(draws, interval, interval.type = c("hpd", "quantile")) {
  interval.type <- match.arg(interval.type)
  draws <- as.numeric(draws)
  draws <- draws[is.finite(draws)]
  if (length(draws) < 2L) {
    return(c(NA_real_, NA_real_))
  }
  if (interval.type == "hpd") {
    as.numeric(coda::HPDinterval(coda::mcmc(draws), prob = interval)[1, ])
  } else {
    alpha <- 1 - interval
    as.numeric(stats::quantile(
      draws,
      probs = c(alpha / 2, 1 - alpha / 2),
      names = FALSE,
      type = 8
    ))
  }
}

# Harmonic mean
hmean <- function(x) {
  return(1 / mean(1 / x))
}

# Compute the LPML criteria:
LPML <- function(loglik) {
  lik <- apply(loglik, 2, exp)
  CPO <- apply(lik, 2, hmean)
  LPML <- sum(log(CPO))
  aLPML <- mean(log(CPO))
  return(c(LPML, aLPML))
}

# Compute the DIC criteria:
DIC <- function(loglik) {
  D <- apply(-2 * loglik, 1, sum)
  pD <- 0.5 * stats::var(D)
  DIC <- mean(D) + pD
  return(matrix(c(DIC, pD), ncol = 2))
}

# Compute the WAIC criteria:
WAIC <- function(loglik) {
  lpd <- sum(log(apply(exp(loglik), 2, mean)))
  pD <- sum(apply(loglik, 2, stats::var))
  WAIC <- lpd - pD
  return(matrix(c(WAIC, pD), ncol = 2))
}

#' Internal: Calculate the posterior mode
#'
#' @keywords internal
.mode <- function(ext) {
  f <- density(ext)
  pmode <- f$x[which.max(f$y)]
  return(pmode)
}

#' Internal: AFT basis evaluator aligned with Stan transformed parameters
#'
#' @param time_aft Numeric vector of AFT residual-scale times y = log(t) - eta.
#' @param tau_a Lower bound used by the fitted model.
#' @param tau_b Upper bound used by the fitted model.
#' @param P Power-basis transformation matrix from \code{pw.basis()}.
#'
#' @return \code{list(g = ..., G = ...)} where \code{g} is hazard basis and
#'   \code{G} is cumulative-hazard basis on the AFT residual scale.
#' @keywords internal
.spbp_aft_basis <- function(time_aft, tau_a, tau_b, P) {
  m <- ncol(P)
  if (m <= 0L) {
    stop("invalid power-basis matrix 'P'", call. = FALSE)
  }
  range_aft <- tau_b - tau_a
  if (!is.finite(range_aft) || range_aft <= 0) {
    stop("invalid AFT range: need tau_b > tau_a", call. = FALSE)
  }

  u <- (as.numeric(time_aft) - tau_a) / range_aft
  n <- length(u)
  b <- matrix(0, nrow = n, ncol = m)
  B <- matrix(0, nrow = n, ncol = m)
  for (j in seq_len(m)) {
    b[, j] <- u^(j - 1)
    B[, j] <- u^j / j
  }

  list(
    g = (b %*% P) / range_aft,
    G = (B %*% P)
  )
}

#' AFT residual-scale range tau_a(beta), tau_b(beta) from training data.
#'
#' @noRd
.spbp_aft_tau_range <- function(beta, X_train, log_time_train, tau_a = NULL, tau_b = NULL) {
  if (!is.null(X_train) && !is.null(log_time_train)) {
    w_train <- log_time_train - as.numeric(X_train %*% beta)
    list(tau_a = min(w_train), tau_b = max(w_train))
  } else {
    list(tau_a = tau_a, tau_b = tau_b)
  }
}

#' Symmetric matrix inverse: Cholesky if SPD, else full-rank QR solve, else pivoted-QR pseudoinverse.
#'
#' @noRd
.spbp_sym_inv <- function(A, tol = 1e-10) {
  if (!is.matrix(A) || nrow(A) == 0L || ncol(A) == 0L || nrow(A) != ncol(A)) {
    stop("symmetric inverse requires a square non-empty matrix", call. = FALSE)
  }
  if (!is.numeric(A) || any(!is.finite(A))) {
    stop("symmetric inverse requires a finite numeric matrix", call. = FALSE)
  }
  A <- (A + t(A)) / 2
  n <- nrow(A)
  Rch <- tryCatch(chol(A), error = function(e) NULL)
  if (!is.null(Rch)) {
    return(chol2inv(Rch))
  }
  qrA <- tryCatch(
    qr(A, tol = tol, LAPACK = TRUE),
    error = function(e) NULL
  )
  if (is.null(qrA)) {
    stop("QR decomposition failed for symmetric inverse", call. = FALSE)
  }
  if (qrA$rank == n) {
    return(qr.solve(qrA, diag(n)))
  }
  .spbp_sym_pinv_piv_qr(A, tol = tol)
}

#' Moore-Penrose inverse for symmetric matrices via pivoted QR.
#'
#' @noRd
.spbp_sym_pinv_piv_qr <- function(A, tol = 1e-10) {
  A <- (A + t(A)) / 2
  n <- nrow(A)
  qrA <- qr(A, tol = tol, LAPACK = TRUE)
  r <- qrA$rank
  if (r == 0L) {
    return(matrix(0, n, n))
  }
  Rr <- qr.R(qrA)[seq_len(r), seq_len(r), drop = FALSE]
  Q <- qr.Q(qrA)[, seq_len(r), drop = FALSE]
  Rinv <- qr.solve(qr(Rr), diag(r))
  mid <- rbind(Rinv, matrix(0, n - r, r))
  res <- matrix(0, n, n)
  res[qrA$pivot[seq_len(r)], ] <- mid %*% t(Q)
  res
}

#' Quadratic form g' V g for symmetric V without explicit inversion when SPD.
#'
#' @noRd
.spbp_sym_quad_form <- function(V, g) {
  g <- as.numeric(g)
  V <- (V + t(V)) / 2
  Rch <- tryCatch(chol(V), error = function(e) NULL)
  if (!is.null(Rch)) {
    z <- as.numeric(Rch %*% g)
    return(sum(z * z))
  }
  as.numeric(t(g) %*% V %*% g)
}

#' Extract the negative Hessian on (beta, gamma), dropping alpha if present.
#'
#' @noRd
.spbp_hessian_neg_beta_gamma <- function(object) {
  q <- length(object$coefficients)
  m <- length(object$bp.param)
  H_raw <- object$hessian
  if (is.null(H_raw) || !is.matrix(H_raw)) {
    return(list(H_neg = NULL, q = q, m = m))
  }

  H_neg <- -as.matrix(H_raw)
  nr <- nrow(H_neg)
  target <- q + m
  full_w_alpha <- 1L + target
  if (nr == full_w_alpha) {
    H_neg <- H_neg[-1L, -1L, drop = FALSE]
  } else if (nr != target) {
    return(list(H_neg = NULL, q = q, m = m))
  }

  list(H_neg = H_neg, q = q, m = m)
}

#' Diagnostics for the gamma (Bernstein polynomial) information block.
#'
#' @param kappa_threshold Condition numbers above this are treated as unstable.
#' @noRd
.spbp_gamma_information_diagnostics <- function(object, kappa_threshold = 1e10) {
  hx <- .spbp_hessian_neg_beta_gamma(object)
  m <- hx$m
  unavailable <- function(reason) {
    list(
      stable = FALSE,
      available = FALSE,
      kappa_gamma = NA_real_,
      rank_gamma = NA_integer_,
      n_gamma = m,
      reason = reason
    )
  }
  if (is.null(hx$H_neg) || m == 0L) {
    return(unavailable("gamma information matrix is unavailable"))
  }

  q <- hx$q
  C <- if (q == 0L) {
    hx$H_neg
  } else {
    hx$H_neg[(q + 1L):(q + m), (q + 1L):(q + m), drop = FALSE]
  }

  if (!is.matrix(C) || !is.numeric(C) || nrow(C) == 0L || ncol(C) == 0L) {
    return(unavailable("gamma information matrix is not a valid numeric matrix"))
  }
  if (nrow(C) != m || ncol(C) != m) {
    return(unavailable("gamma information matrix has invalid dimensions"))
  }
  if (any(!is.finite(C))) {
    return(unavailable("gamma information matrix contains non-finite values"))
  }

  qrC <- tryCatch(qr(C), error = function(e) NULL)
  if (is.null(qrC)) {
    return(unavailable("QR decomposition of the gamma information matrix failed"))
  }
  rank_gamma <- qrC$rank
  kappa_gamma <- tryCatch(kappa(C), error = function(e) NA_real_)
  stable <- is.finite(kappa_gamma) &&
    kappa_gamma <= kappa_threshold &&
    !is.na(rank_gamma) &&
    rank_gamma >= m

  list(
    stable = stable,
    available = TRUE,
    kappa_gamma = kappa_gamma,
    rank_gamma = rank_gamma,
    n_gamma = m,
    reason = if (stable) {
      NA_character_
    } else {
      "gamma information matrix is ill-conditioned"
    }
  )
}

#' Warn when survfit delta-method bands may be unreliable due to unstable gamma information.
#'
#' @noRd
.spbp_warn_gamma_unstable_survfit <- function(object, kappa_gamma, reason = NULL) {
  m <- length(object$bp.param)
  if (!is.null(reason) && !is.na(reason) && identical(reason, "gamma information matrix contains non-finite values")) {
    warning(
      "Bernstein-polynomial (gamma) information matrix contains non-finite values ",
      "(degree = ", m, "); ",
      "delta-method survival standard errors and confidence bands are unavailable. ",
      "Try refitting with a lower Bernstein degree ",
      "(argument `degree` to bpph(), bppo(), bpaft(), or spbp()).",
      call. = FALSE
    )
    return(invisible(NULL))
  }
  kappa_txt <- if (is.finite(kappa_gamma)) signif(kappa_gamma, 3) else "NA"
  warning(
    "Bernstein-polynomial (gamma) information matrix is ill-conditioned ",
    "(kappa = ", kappa_txt, ", degree = ", m, "); ",
    "delta-method survival standard errors and confidence bands may be unreliable. ",
    "Try refitting with a lower Bernstein degree ",
    "(argument `degree` to bpph(), bppo(), bpaft(), or spbp()).",
    call. = FALSE
  )
}

#' Internal: compute survival CI bands (survfit-style)
#'
#' @keywords internal
#' @importFrom stats qnorm
.survfit_confint <- function(p, se, logse = TRUE, conf.type, conf.int, selow, ulimit = TRUE) {
  zval <- qnorm(1 - (1 - conf.int) / 2, 0, 1)
  if (missing(selow)) {
    scale <- 1
  } else {
    scale <- ifelse(selow == 0, 1, selow / se)
  }
  if (!logse) {
    se <- ifelse(se == 0, 0, se / p)
  }
  if (conf.type == "plain") {
    se2 <- se * p * zval
    if (ulimit) {
      list(lower = pmax(p - se2 * scale, 0), upper = pmin(p +
        se2, 1))
    } else {
      list(lower = pmax(p - se2 * scale, 0), upper = p +
        se2)
    }
  } else if (conf.type == "log") {
    xx <- ifelse(p == 0, NA, p)
    se2 <- zval * se
    temp1 <- exp(log(xx) - se2 * scale)
    temp2 <- exp(log(xx) + se2)
    if (ulimit) {
      list(lower = temp1, upper = pmin(temp2, 1))
    } else {
      list(lower = temp1, upper = temp2)
    }
  } else if (conf.type == "log-log") {
    xx <- ifelse(p == 0 | p == 1, NA, p)
    se2 <- zval * se / log(xx)
    temp1 <- exp(-exp(log(-log(xx)) - se2 * scale))
    temp2 <- exp(-exp(log(-log(xx)) + se2))
    list(lower = temp1, upper = temp2)
  } else if (conf.type == "logit") {
    xx <- ifelse(p == 0, NA, p)
    se2 <- zval * se * (1 + xx / (1 - xx))
    temp1 <- 1 - 1 / (1 + exp(log(p / (1 - p)) - se2 * scale))
    temp2 <- 1 - 1 / (1 + exp(log(p / (1 - p)) + se2))
    list(lower = temp1, upper = temp2)
  } else if (conf.type == "arcsin") {
    xx <- ifelse(p == 0, NA, p)
    se2 <- 0.5 * zval * se * sqrt(xx / (1 - xx))
    list(
      lower = (sin(pmax(0, asin(sqrt(xx)) - se2 * scale)))^2,
      upper = (sin(pmin(pi / 2, asin(sqrt(xx)) + se2)))^2
    )
  } else {
    stop("invalid conf.int type")
  }
}

#' @keywords internal
#' @noRd
.spbp_summary_message <- function(x) {
  cls <- class(x)[1]
  is_bayes <- grepl("bayes$", cls)
  prefix <- if (is_bayes) "Bayesian " else ""
  model <- sub("^summary\\.bp", "", sub("\\.(mle|bayes)$", "", cls))
  model_label <- switch(model,
    ph = "PH",
    po = "PO",
    aft = "AFT",
    toupper(model)
  )
  paste0(prefix, "Bernstein ", model_label, " model: \n")
}

#' @keywords internal
#' @noRd
.spbp_drop_all_na_cols <- function(x) {
  if (is.null(x)) {
    return(x)
  }
  if (is.matrix(x)) {
    if (ncol(x) == 0L) {
      return(x)
    }
    keep <- apply(x, 2, function(col) !all(is.na(col)))
    return(x[, keep, drop = FALSE])
  }
  if (!is.data.frame(x) || ncol(x) == 0L) {
    return(x)
  }
  keep <- vapply(x, function(col) !all(is.na(col)), logical(1))
  x[, keep, drop = FALSE]
}

#' @keywords internal
#' @noRd
.spbp_interval_colnames <- function(interval_colnames) {
  lower_col <- grep("^lower", interval_colnames, value = TRUE)[1]
  pct <- as.numeric(sub("^lower \\.", "", sub("HPD$", "", lower_col)))
  level <- pct / 100
  alpha <- (1 - level) / 2
  c(
    paste0(round(100 * alpha, 1), "%"),
    paste0(round(100 * (1 - alpha), 1), "%")
  )
}
