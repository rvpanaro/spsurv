#' Model.matrix method for fitted spbp models
#'
#' @export
#' @param object an object of class `spbp`, see \code{\link[spsurv]{spbp}}.
#' @method model.matrix spbp
#' @description Model.matrix of a fitted \code{\link[spsurv]{spbp}} model.
#' @param data data.frame object.
#' @param ... arguments inherent from \code{\link[stats]{model.matrix}}.
#' @return The explanatory variables matrix.
#' @seealso \code{\link[spsurv]{spbp}}, \code{\link[stats]{model.matrix}}
#' @importFrom stats model.matrix
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
  function(x, ...) {
    x$features
  }
