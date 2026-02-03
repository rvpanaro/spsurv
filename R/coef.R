#---------------------------------------------
#' Estimated regression coefficients
#'
#' @aliases coef.spbp
#' @export
#' @param object an object of the class spbp
#' @param summary posterior summary if method ="bayes" in x
#' @param ... further arguments passed to or from other methods
#' @return  the estimated regression coefficients
#' @importFrom stats coef
#'
#'
coef.spbp <- function(object, summary = c("mean", "median", "mode"), ...) {
  if (object$call$approach == "mle") {
    return(object$coefficients)
  } else if (!is.null(object$posterior$beta)) {
    summary_flag <- match.arg(summary)

    if (summary_flag == "mean") {
      return(apply(object$posterior$beta, 2, mean))
    } else if (summary_flag == "median") {
      return(apply(object$posterior$beta, 3, median))
    } else {
      return(apply(object$posterior$beta, 3, mode))
    }
  } else {
    return(NULL)
  }
}
