#' @export
#' @method survfit spbp
#' @title BP-based model survival curves
#' @description Compute survival curves for a fitted \code{\link[spsurv:spbp]{spbp}} model.
#'
#' @param formula An object of class \code{"spbp"} returned by \code{\link[spsurv:spbp]{spbp}}.
#' @param newdata Optional data frame used to obtain survival curves for specific covariate values.
#' @param times Optional numeric vector of time points at which to return estimates.
#' @param se.fit Logical; if \code{TRUE}, compute standard errors.
#' @param interval Confidence level for intervals (e.g. \code{0.95}).
#' @param type Character; confidence interval transformation. One of \code{"log"}, \code{"log-log"}, or \code{"plain"}.
#' @param ... Further arguments (currently ignored or reserved for future use).
#' @importFrom coda mcmc
#' @importFrom stats model.matrix
#' @return An object of class \code{"survfit"}.
#' @seealso \code{\link[spsurv:spbp]{spbp}}, \code{\link[survival:survfit]{survfit}}.
#' @examples
#' library(spsurv)
#' data(veteran, package = "survival")
#' fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran)
#' survfit(fit)
#'
survfit.spbp <- function(formula, newdata, times,
                         se.fit = TRUE, interval = .95,
                         type = c("log", "log-log", "plain"),
                         ...) {
  x <- formula
  type <- match.arg(type)

  data <- eval(x$call$data, envir = environment(x$terms))

  coxfit <- coxph(x$formula, data = data, model = TRUE)
  p <- length(x$coefficients)

  if (missing(newdata)) {
    km <- survival::survfit(formula = coxfit, ...)
    X <- t(x$means)
  } else {
    km <- survival::survfit(formula = coxfit, newdata = newdata, ...)
    X <- matrix(model.matrix(object = x$formula[-2], xlev = x$xlevels, data = newdata)[, -1], ncol = p)
  }

  if (!missing(times)) {
    if (!inherits(times, "Surv")) {
      stop("times is not a 'Surv' object")
    } else {
      km$time <- sort(times[, 1], decreasing = F)
      km$n.event <- times[, 2][order(times[, 1], decreasing = F)]
      km$n.censor <- 1 - km$n.event
      km$n.censor <- 1 - km$n.event
      km$n.risk <- length(times[, 1]) - cumsum(km$n.event) + 1
    }
  }

  m <- length(x$bp.param)
  cumhaz <- NULL

  if (x$call$approach == "mle") {
    if (p == 0) {
      beta <- 0
    } else {
      beta <- x$coefficients
    }
    gamma <- x$bp.param

    exp_eta <- exp(X %*% beta) %>% as.vector()
    var <- vcov.spbp(formula, bp.param = TRUE)

    if (x$call$model == "ph") {
      G <- matrix(sapply(1:m, function(k) pbeta(km$time / x$tau_b, k, m - k + 1)), nrow = length(km$time))

      cumhaz <- (G %*% gamma %>% as.vector()) %o% exp_eta
      surv <- exp(-cumhaz)
      std.err <- cumhaz

      for (i in 1:nrow(cumhaz)) {
        for (j in 1:ncol(cumhaz)) {
          dbeta <- (cumhaz[i, j] * X[j, ])
          dgamma <- (gamma * G[i, ] * exp_eta[j])

          if (p > 0) {
            grad <- as.matrix(c(dbeta, dgamma))
          } else {
            grad <- as.matrix(c(dgamma))
          }

          std.err[i, j] <-
            as.numeric(t(grad) %*% var %*% grad)
        }
      }

      std.err <- sqrt(std.err)
    } else if (x$call$model == "po") {
      G <- matrix(sapply(1:m, function(k) pbeta(km$time / x$tau_b, k, m - k + 1)), nrow = length(km$time))

      odds <- (G %*% gamma %>% as.vector()) %o% exp_eta
      cumhaz <- log(1 + odds)
      surv <- exp(-cumhaz)
      std.err <- odds

      for (i in 1:nrow(odds)) {
        for (j in 1:ncol(odds)) {
          dbeta <- 1 / (1 + odds[i, j]) * (odds[i, j] * X[j, ])
          dgamma <- 1 / (1 + odds[i, j]) * (gamma * G[i, ] * exp_eta[j])

          if (p > 0) {
            grad <- as.matrix(c(dbeta, dgamma))
          } else {
            grad <- as.matrix(c(dgamma))
          }

          std.err[i, j] <-
            as.numeric(t(grad) %*% var %*% grad)
        }
      }

      std.err <- sqrt(std.err)
    } else {
      time <- log(km$time %o% (1 / exp_eta))

      cumhaz <- matrix(ncol = nrow(X), nrow = length(km$time))
      haz <- cumhaz
      g <- list()
      G <- list()

      for (j in 1:nrow(X)) {
        g[[j]] <- matrix(sapply(1:m, function(k) dbeta((time[, j] - x$tau_a) / (x$tau_b - x$tau_a), k, m - k + 1)) / (x$tau_b - x$tau_a), nrow = length(km$time))

        haz[, j] <- (g[[j]] %*% gamma) %>%
          as.vector()

        G[[j]] <- matrix(sapply(1:m, function(k) pbeta((time[, j] - x$tau_a) / (x$tau_b - x$tau_a), k, m - k + 1)), nrow = length(km$time))

        cumhaz[, j] <- (G[[j]] %*% gamma) %>%
          as.vector()
      }

      surv <- exp(-cumhaz)
      std.err <- cumhaz

      for (i in 1:nrow(cumhaz)) {
        for (j in 1:ncol(cumhaz)) {
          dbeta <- -(haz[i, j] * X[j, ])
          dgamma <- (gamma * G[[j]][i, ])

          if (p > 0) {
            grad <- as.matrix(c(dbeta, dgamma))
          } else {
            grad <- as.matrix(c(dgamma))
          }

          std.err[i, j] <-
            as.numeric(t(grad) %*% var %*% grad)
        }
      }
      std.err <- sqrt(std.err)
    }

    km$surv <- surv[, 1:nrow(X)]
    km$cumhaz <- cumhaz[, 1:nrow(X)]
    km$std.err <- std.err[, 1:nrow(X)]

    ci <-
      .survfit_confint(
        p = surv,
        se = std.err,
        logse = km$logse,
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
    cumhaz <- array(dim = c(length(km$time), nrow(X), n.samp))
    G <- sapply(1:m, function(k) pbeta(km$time / x$tau_b, k, m - k + 1))

    if (x$call$model == "ph") {
      for (i in 1:n.samp) {
        exp_eta[i, ] <- exp(X %*% beta[i, ]) %>% as.vector()
        cumhaz[, , i] <- (G %*% gamma[i, ] %>% as.vector()) %o% exp_eta[i, ]
      }
    } else if (x$call$model == "po") {
      odds <- cumhaz
      for (i in 1:n.samp) {
        exp_eta[i, ] <- exp(X %*% beta[i, ]) %>% as.vector()
        odds[, , i] <- (G %*% gamma[i, ] %>% as.vector()) %o% exp_eta[i, ]
      }
      cumhaz <- log(1 + odds)
    } else {
      time <- array(dim = c(length(km$time), nrow(X), n.samp))

      for (i in 1:n.samp) {
        exp_eta[i, ] <- exp(X %*% beta[i, ]) %>% as.vector()
        time[, , i] <- log(km$time %o% (1 / exp_eta[i, ]))
      }

      G <- array(dim = c(length(km$time), m, nrow(X), n.samp))
      for (j in 1:n.samp) {
        for (i in 1:nrow(X)) {
          G[, , i, j] <- sapply(1:m, function(k) pbeta((time[, i, j] - x$tau_a) / (x$tau_b - x$tau_a), k, m - k + 1))

          cumhaz[, i, j] <- (G[, , i, j] %*% gamma[j, ]) %>%
            as.vector()
        }
      }
    }

    surv <- apply(exp(-cumhaz), c(1, 2), mean)
    std.err <- apply(exp(-cumhaz), c(1, 2), sd)

    km$surv <- surv[, 1:nrow(X)]
    km$std.err <- std.err[, 1:nrow(X)]

    ci <- NULL
    if (type == "plain") {
      ci$lower <- apply(exp(-cumhaz), c(1, 2), function(i) {
        coda::HPDinterval(mcmc(i))[, 1]
      })
      ci$upper <- apply(exp(-cumhaz), c(1, 2), function(i) {
        coda::HPDinterval(mcmc(i))[, 2]
      })
    } else if (type == "log") {
      ci$lower <- apply(-cumhaz, c(1, 2), function(i) {
        coda::HPDinterval(mcmc(i))[, 1]
      }) %>% exp()
      ci$upper <- apply(-cumhaz, c(1, 2), function(i) {
        coda::HPDinterval(mcmc(i))[, 2]
      }) %>% exp()
    } else if (type == "log-log") {
      ci$lower <- -exp(apply(log(cumhaz), c(1, 2), function(i) {
        coda::HPDinterval(mcmc(i))[, 2]
      })) %>% exp()
      ci$upper <- -exp(apply(log(cumhaz), c(1, 2), function(i) {
        coda::HPDinterval(mcmc(i))[, 1]
      })) %>% exp()
    }

    cumhaz <- apply(cumhaz, c(1, 2), mean)
    km$cumhaz <- cumhaz[, 1:nrow(X)]
  }

  km$lower <- ci$lower[, 1:nrow(X)]
  km$upper <- ci$upper[, 1:nrow(X)]
  km$conf.type <- type
  km$conf.int <- interval
  km$call <- match.call() ## Call

  if (nrow(X) > 1) {
    colnames(km$surv) <- 1:nrow(X)
    colnames(km$cumhaz) <- 1:nrow(X)
    colnames(km$std.err) <- 1:nrow(X)
    colnames(km$lower) <- 1:nrow(X)
    colnames(km$upper) <- 1:nrow(X)
  }

  class(km) <- c("survfitbp", "survfitcox", "survfit")

  return(km)
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
#' data("veteran")
#'
#' fit <- bpph(Surv(time, status) ~ karno + factor(celltype),
#'   data = veteran
#' )
#'
#' residuals(fit)
residuals.spbp <- function(object, type = c("martingale", "deviance", "coobject-snell"), ...) {
  type <- match.arg(type)

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
    exp_eta <- exp(X %*% beta) %>% as.vector()

    if (object$call$model == "ph") {
      G <- matrix(sapply(1:m, function(k) pbeta(object$y[, 1] / object$tau_b, k, m - k + 1)), nrow = length(object$y[, 1]))
      cumhaz <- (G %*% gamma %>% as.vector()) * exp_eta
    } else if (object$call$model == "po") {
      G <- matrix(sapply(1:m, function(k) pbeta(object$y[, 1] / object$tau_b, k, m - k + 1)), nrow = length(object$y[, 1]))
      odds <- (G %*% gamma %>% as.vector()) * exp_eta
      cumhaz <- log(1 + odds)
    } else {
      time <- log(object$y[, 1]) - (object %*% beta)
      cumhaz <- matrix(ncol = nrow(object), nrow = length(object$y[, 1]))
      G <- matrix(nrow = length(object$y[, 1]), ncol = m)

      for (k in 1:m) {
        G[, k] <- pbeta((time - object$tau_a) / (object$tau_b - object$tau_a), k, m - k + 1)
      }
      cumhaz <- G %*% gamma %>% as.vector()
    }
  } else {
    beta <- object$posterior$beta
    gamma <- object$posterior$gamma
    n.samp <- nrow(beta)
    exp_eta <- matrix(nrow = n.samp, ncol = length(object$y[, 1]))
    cumhaz <- array(dim = c(n.samp, length(object$y[, 1])))
    G <- sapply(1:m, function(k) pbeta(object$y[, 1] / object$tau_b, k, m - k + 1))

    if (object$call$model == "ph") {
      for (i in 1:n.samp) {
        exp_eta[i, ] <- exp(X %*% beta[i, ]) %>% as.vector()
        cumhaz[i, ] <- (G %*% gamma[i, ] %>% as.vector()) * exp_eta[i, ]
      }
    } else if (object$call$model == "po") {
      odds <- cumhaz
      for (i in 1:n.samp) {
        exp_eta[i, ] <- exp(X %*% beta[i, ]) %>% as.vector()
        odds[i, ] <- (G %*% gamma[i, ] %>% as.vector()) * exp_eta[i, ]
      }
      cumhaz <- log(1 + odds)
    } else {
      time <- array(dim = c(n.samp, length(object$y[, 1])))

      for (i in 1:n.samp) {
        exp_eta[i, ] <- exp(X %*% beta[i, ]) %>% as.vector()
        time[i, ] <- log(object$y[, 1] / exp_eta[i, ])
      }

      G <- array(dim = c(n.samp, length(object$y[, 1]), m))
      for (i in 1:n.samp) {
        G[i, , ] <- sapply(1:m, function(k) pbeta((time[i, ] - object$tau_a) / (object$tau_b - object$tau_a), k, m - k + 1))

        cumhaz[i, ] <- (G[i, , ] %*% gamma[i, ]) %>% as.vector()
      }
    }
    cumhaz <- apply(cumhaz, 2, mean)
  }

  cumhaz <- cumhaz %>% as.vector()
  names(cumhaz) <- names(object$y[, 1])
  delta <- object$y[, 2]

  if (type == "coobject-snell") {
    return(cumhaz)
  } else if (type == "martingale") {
    return(delta - cumhaz)
  } else {
    m <- delta - cumhaz
    sign(m) * sqrt(-2 * (m + delta * log(delta - m)))
  }
}
