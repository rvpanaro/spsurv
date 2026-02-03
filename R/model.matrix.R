#' Model.matrix method for fitted spbp models
#'
#' @export
#' @param object an object of class `spbp`, see \code{\link[spsurv]{spbp}}.
#' @method model.matrix spbp
#' @description Model.matrix of a fitted \code{\link[spsurv]{spbp}} model.
#' @return The model matrix.
#' @seealso \code{\link[spsurv]{spbp}}, \code{\link[stats]{model.matrix}}
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
#' model.matrix(fit)
model.matrix.spbp <-
  function(object, ...) {
    object$features
  }
