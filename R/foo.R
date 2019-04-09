#' Bernstein Polynomials for proportional odds model (base odds estimation)
#'
#' @export
#' @param formula a formula object, with time-to-event data on the left side of ~ and explanatory terms on the right.
#' @param m Bernstein Polynomial degree
#' @param tau must be greater than times maximum value observed, this value is used for correction purposes.
#' @param optimize bayes or mle inference method, by default optimize = FALSE.
#' @param data a data.frame with variables named in the formula.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'

bpph <- function(...){
  spbpArgs <- list(...)

  if (length(stanArgs)) {
    formals <- names(formals(spbp::spbp))
    aux <- pmatch(names(spbpArgs), formals, nomatch = 0)

    if (any(aux == 0))
      stop(gettextf("Argument %s not matched", names(spbpArgs)[aux==0]))
  }
  else{
    stop("Arguments not matched")
  }
  spbp(...)
}

#' Bernstein Polynomials for proportional hazards model (base risk estimation)
#'
#' @export
#' @param formula a formula object, with time-to-event data on the left side of ~ and explanatory terms on the right.
#' @param m Bernstein Polynomial degree
#' @param tau must be greater than times maximum value observed, this value is used for correction purposes.
#' @param optimize bayes or mle inference method, by default optimize = FALSE.
#' @param data a data.frame with variables named in the formula.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'

bpph <- function(...){
  spbpArgs <- list(...)

  if (length(stanArgs)) {
    formals <- list(names(formals(spbp::spbp)),
                     names(formals(rstan::stan)),
                     names(formals(rstan::optimizing)))

    aux <- pmatch(names(spbpArgs), formals, nomatch = 0)

    if (any(aux == 0))
      stop(gettextf("Argument %s not matched", names(spbpArgs)[aux==0]))
  }
  else{
    stop("Arguments not matched")
  }
  spbp(...)
}
