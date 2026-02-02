#'
#' @export
terms.inner <- function(x) {
  if (class(x) == "formula") {
    c(terms.inner(x[[2]]), terms.inner(x[[3]]))
  } else if (class(x) == "call" && (x[[1]] != as.name("$") &&
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
  pD <- 0.5 * var(D)
  DIC <- mean(D) + pD
  return(matrix(c(DIC, pD), ncol = 2))
}

# Compute the WAIC criteria:
WAIC <- function(loglik) {
  lpd <- sum(log(apply(exp(loglik), 2, mean)))
  pD <- sum(apply(loglik, 2, var))
  WAIC <- lpd - pD
  return(matrix(c(WAIC, pD), ncol = 2))
}

## Calculate the posterior mode
mode <- function(ext) {
  f <- density(ext)
  pmode <- f$x[which.max(f$y)]
  return(pmode)
}
