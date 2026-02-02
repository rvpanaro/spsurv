#' @export
#' @method survfit spbp
#' @title BP based models survfit.
#' @description Survival for a fitted \code{\link[spsurv]{spbp}} model.
#' @param x an object of class `spbp` result of a \code{\link[spsurv]{spbp}} fit.
#' @param ... further arguments passed to survfit.
#' @seealso \code{\link[spsurv]{spbp}}, \code{\link[survfit]{spbp}}.
#' @examples
#'
#' library("spsurv")
#' data("veteran")
#'
#' fit <- bpph(Surv(time, status) ~ karno + factor(celltype),
#'   data = veteran
#' )
#'
#' survfit(fit)
survfit.spbp <- function(x, newdata, times,
                         se.fit = TRUE, interval = .95,
                         type = c("log", "log-log", "plain"),
                         ...) {
  formula <- x$formula
  suppressMessages(eval(expr = x$call$data, envir = environment(x$terms)) %>%
    attach())
  type <- match.arg(type)

  coxfit <- coxph(formula)
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

    exp_eta <-
      exp(X %*% beta) %>%
      as.vector()

    var <- vcov.spbp(x, bp.param = TRUE)

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
      survival:::survfit_confint(
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
    n.samp <- nrow(beta)
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
#' @param x an object of class `spbp` result of a \code{\link[spsurv]{spbp}} fit.
#' @param type type of residuals, default is "cox-snell"
#' @param ... further arguments passed to or from other methods
#' @seealso \code{\link[spsurv]{spbp}}, \code{\link[survfit]{spbp}}.
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
residuals.spbp <- function(x, type = c("martingale", "deviance", "cox-snell")) {
  type <- match.arg(type)

  p <- length(x$coefficients)

  if (p > 0) {
    X <- model.matrix(x)
  } else {
    X <- t(x$means)
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

    if (x$call$model == "ph") {
      G <- matrix(sapply(1:m, function(k) pbeta(x$y[, 1] / x$tau_b, k, m - k + 1)), nrow = length(x$y[, 1]))
      cumhaz <- (G %*% gamma %>% as.vector()) * exp_eta
    } else if (x$call$model == "po") {
      G <- matrix(sapply(1:m, function(k) pbeta(x$y[, 1] / x$tau_b, k, m - k + 1)), nrow = length(x$y[, 1]))
      odds <- (G %*% gamma %>% as.vector()) * exp_eta
      cumhaz <- log(1 + odds)
    } else {
      time <- log(x$y[, 1]) - (X %*% beta)
      cumhaz <- matrix(ncol = nrow(X), nrow = length(x$y[, 1]))
      G <- matrix(nrow = length(x$y[, 1]), ncol = m)

      for (k in 1:m) {
        G[, k] <- pbeta((time - x$tau_a) / (x$tau_b - x$tau_a), k, m - k + 1)
      }
      cumhaz <- G %*% gamma %>% as.vector()
    }
  } else {
    beta <- x$posterior$beta
    gamma <- x$posterior$gamma
    n.samp <- nrow(beta)
    exp_eta <- matrix(nrow = n.samp, ncol = length(x$y[, 1]))
    cumhaz <- array(dim = c(n.samp, length(x$y[, 1])))
    G <- sapply(1:m, function(k) pbeta(x$y[, 1] / x$tau_b, k, m - k + 1))

    if (x$call$model == "ph") {
      for (i in 1:n.samp) {
        exp_eta[i, ] <- exp(X %*% beta[i, ]) %>% as.vector()
        cumhaz[i, ] <- (G %*% gamma[i, ] %>% as.vector()) * exp_eta[i, ]
      }
    } else if (x$call$model == "po") {
      odds <- cumhaz
      for (i in 1:n.samp) {
        exp_eta[i, ] <- exp(X %*% beta[i, ]) %>% as.vector()
        odds[i, ] <- (G %*% gamma[i, ] %>% as.vector()) * exp_eta[i, ]
      }
      cumhaz <- log(1 + odds)
    } else {
      time <- array(dim = c(n.samp, length(x$y[, 1])))

      for (i in 1:n.samp) {
        exp_eta[i, ] <- exp(X %*% beta[i, ]) %>% as.vector()
        time[i, ] <- log(x$y[, 1] / exp_eta[i, ])
      }

      G <- array(dim = c(n.samp, length(x$y[, 1]), m))
      for (i in 1:n.samp) {
        G[i, , ] <- sapply(1:m, function(k) pbeta((time[i, ] - x$tau_a) / (x$tau_b - x$tau_a), k, m - k + 1))

        cumhaz[i, ] <- (G[i, , ] %*% gamma[i, ]) %>% as.vector()
      }
    }
    cumhaz <- apply(cumhaz, 2, mean)
  }

  cumhaz <- cumhaz %>% as.vector()
  names(cumhaz) <- names(x$y[, 1])
  delta <- x$y[, 2]

  if (type == "cox-snell") {
    return(cumhaz)
  } else if (type == "martingale") {
    return(delta - cumhaz)
  } else {
    m <- delta - cumhaz
    sign(m) * sqrt(-2 * (m + delta * log(delta - m)))
  }
}


#' Plot the residuals for the Cox model
#'
#' @param spbp A coxfit object
#' @param type one of "cox-snell", "martingale" or "deviance"
#' @return The residuals of selected type
#' @examples
#' library(survverse)
#' coxph.fit2 <- coxph(Surv(futime, fustat) ~ age + ecog.ps, data = ovarian)
#' ggcoxdiagnostics2(coxph.fit2, type = "deviance")
#' ggcoxdiagnostics2(coxph.fit2, type = "schoenfeld")
ggdiagnostics2 <- function(coxfit, type) {
  type %<>% match.arg(arg = ., choices = c(
    "cox-snell", "martingale",
    "deviance", "schoenfeld"
  ))

  if (!inherits(x = coxfit, what = "coxph")) {
    stop("survival's coxph object required")
  }
  delta <- coxfit$y[, 2]

  if (type == "cox-snell") {
    res <- delta - (coxfit %>% resid())
    KM <- survfit(Surv(res, delta) ~ 1, type = "kaplan-meier")

    ggplot() +
      geom_step(mapping = aes(
        x = KM$time,
        y = -log(KM$surv)
      )) +
      scale_y_continuous(limits = c(0, 4)) +
      xlab("Residuals") +
      ylab("Cumulative risk") +
      ggtitle("Cox model: Cox-Snell residuals") +
      geom_abline(col = c("blue")) +
      theme(text = element_text(size = 12))
  } else if (type %in% c("martingale", "deviance")) {
    if (type == "martingale") {
      res <- (coxfit %>% resid())
    } else {
      res <- (coxfit %>% resid(., type = "deviance"))
    }

    vars <- coxfit %>%
      model.frame() %>%
      select(names(attr(coxfit$terms, "dataClasses"))[attr(coxfit$terms, "dataClasses") == "numeric"])

    for (j in 1:ncol(vars)) {
      print(ggplot() +
        geom_point(aes(
          x = vars[, j], res,
          color = factor(delta)
        )) +
        geom_smooth(aes(x = vars[, j], res),
          method = "loess", se = FALSE, span = 1
        ) +
        scale_color_manual(
          labels = c("Censored", "Non-censored"),
          values = c("blue", "red")
        ) +
        labs(
          title = paste("Functional aspect:", colnames(vars)[j]),
          x = colnames(vars)[j], y = "Martingale residuals", color = "Event"
        ) +
        theme(
          axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          plot.title = element_text(size = 16, face = "bold"),
          legend.position = c(0.9, 0.5)
        ))
    }
  } else if (type == "schoenfeld") {
    zph <- cox.zph(coxfit, transform = "identity")
    plts <- list()

    for (j in 1:ncol(zph$y)) {
      print(ggplot() +
        geom_point(aes(x = zph$x, zph$y[, j])) +
        geom_smooth(aes(x = zph$x, zph$y[, j]),
          method = "loess", se = TRUE, span = 1
        ) +
        labs(
          title = paste0("Schoenfeld residuals", ", p-value ", zph$table[j, "p"] %>% round(2)),
          x = "Time",
          y = paste("Beta(t) for", colnames(zph$y)[j])
        ) +
        theme(
          axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          plot.title = element_text(size = 16, face = "bold"),
          legend.position = c(0.9, 0.5)
        ) +
        geom_hline(yintercept = mean(zph$y[, j]), col = "red", lty = 2))
    }
  }
}
