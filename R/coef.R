#' Estimated regression coefficients
#'
#' @aliases coef.spbp
#' @export
#' @param object an object of the class spbp
#' @param summary posterior summary if method = "bayes" in x
#' @param ... further arguments passed to or from other methods
#' @return  the estimated regression coefficients
#' @importFrom stats coef
coef.spbp <- function(object, summary = c("mean", "median", "mode"), ...) {
  if (identical(object$call$approach, "mle")) {
    return(object$coefficients)
  }
  if (is.null(object$posterior$beta)) {
    return(NULL)
  }

  summary_flag <- match.arg(summary)
  beta_draws <- object$posterior$beta
  switch(summary_flag,
    mean = apply(beta_draws, 2, mean),
    median = apply(beta_draws, 2, median),
    mode = apply(beta_draws, 2, .mode)
  )
}
