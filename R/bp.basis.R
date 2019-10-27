#' Bernstein Polynomials basis calculations
#' @export
#' @param time a vector of times to events.
#' @param degree Bernstein Polynomial degree
#' @param tau must be greater than times maximum value observed, this value is used for correction purposes.
#' @param data a data.frame with variables named in the formula.
#' @return A list containing matrices b and B corresponding BP basis and corresponding tau value used to compute them.

bp.basis <- function(time, degree,  tau = max(time)){
  n <- length(time)

  ## error handling
  if(sum(time >= 0) != n)
    stop("time must be a positive vector.")

  if(degree < 0)
    stop("polynomial degree must be positive.")

  if(!(degree %% 1 == 0))
    stop("polynomial degree must be integer.")

  if(tau <  max(time))
    stop("tau must be greater than the last time.")

  k <- 1:degree
  b <- matrix(NA, n, degree )
  B <- matrix(NA, n, degree )
  y <- time/tau

  b <- sapply(k, function(k){dbeta(y, k, degree - k + 1) / tau})
  B <- sapply(k, function(k) pbeta(y, k, degree - k + 1) )

  # Equivalent to
  # for (i in 1:n){
  #   for(k in 1:degree){
  #     b[i,k] <- dbeta(y[i], k, degree - k + 1) / tau
  #     B[i,k] <- pbeta(y[i], k, degree - k + 1)
  #   }
  # }
  return(list(b = b, B = B, degree = degree, tau = tau))
}

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
####

blockSolve <- function(M, q){

  if(!is.matrix(M))
    stop("M is not a matrix")

  if(!is.numeric(q) || length(q) > 1 || q%%1 != 0)
    stop("q must be an integer")

  if(nrow(M) != ncol(M))
    stop("non square matrix")
  if(q == 0) return(M)

  n <- nrow(M)
  r <- (q+1)

  A <- M[1:q, 1:q];
  B <- M[1:q, r:n];
  C <- M[r:n,1:q];
  D <- M[r:n, r:n];

  S <- matrix(NA, nrow = nrow(M), ncol = ncol(M))

  ## aux op
  invA <- solve(A)
  invPsis <-  MASS::ginv(D - C %*% invA %*% B)

  ##S11
  S[1:q, 1:q] <- invA + invA %*% B %*% invPsis %*% C %*% invA
  ##S12
  S[1:q, r:n] <- (-invA) %*% B %*% invPsis
  ##S21
  S[r:n,1:q] <- (-invPsis) %*% C %*% invA
  ##S22
  S[r:n, r:n] <- invPsis

  return(S)
}

solveAny <- function(A){
  tol <- .Machine$double.eps ## solve default tolerance
  class(S) <- "try-error"
  while(class(S) == "try-error"){
    S <- try(qr.solve(A, tol = tol), silent = T)
    tol <- tol^(1.00001)
  }
  return(S)
}

## Singular Value Decompositition
solveSVD <- function(A){
  svd <- svd(A)
  v <- svd$v
  u <- svd$u
  zeros <- svd$d == 0
  pseudod <- svd$d
  pseudod[!zeros] <- svd$d[!zeros] <-  1/svd$d[!zeros]
  S <- v %*% diag(pseudod) %*% t(u)
  return(S)
}

read_prior <- function(prior){
  aux <- unlist(strsplit(prior, "\\("))
  dist <- aux[1]
  aux2 <- unlist(strsplit(aux[2], "\\)"))[1]
  val <- unlist(strsplit(aux2, "\\,"))
  return(c(dist, val))
}

