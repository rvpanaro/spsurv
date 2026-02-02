#---------------------------------------------
#' Generic S3 method coef
#' @aliases coef
#' @export
#' @param x a fitted model object
#' @param ... further arguments passed to or from other methods.
#' @return the estimated regression coefficients
#'
coef <- function(x, ...) UseMethod("coef")

#---------------------------------------------
#' Estimated regression coefficients
#'
#' @aliases coef.spbp
#' @rdname coef-methods
#' @method coef spbp
#' @export
#' @export coef
#' @param x an object of the class spbp
#' @param ... further arguments passed to or from other methods
#' @return  the estimated regression coefficients
#'
#'
coef.spbp <- function(x, summary = c("mean", "median", "mode"), ...) {
  if (x$call$approach == "mle") {
    return(x$coefficients)
  } else if (!is.null(x$posterior$beta)) {
    summary_flag <- match.arg(summary)

    if (summary_flag == "mean") {
      return(apply(x$posterior$beta, 2, mean))
    } else if (summary_flag == "median") {
      return(apply(x$posterior$beta, 3, median))
    } else {
      return(apply(x$posterior$beta, 3, spsurv:::mode))
    }
  } else {
    return(NULL)
  }
}
