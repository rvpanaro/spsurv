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
  val <- unlist(strsplit(aux2, "\\,"))
  return(c(dist, val))
}

# Harmonic mean
hmean <- function(x) {
  return(1 / mean(1 / x))
}

#  Compute the LPML criteria:
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
