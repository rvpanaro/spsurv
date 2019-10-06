#' Model.matrix method for fitted spbp models
#'
#' @export
#' @method model.matrix spbp
#' @description Model.matrix of a fitted \code{\link[spsurv]{spbp}} model.
#' @param object the result of a \code{\link[spsurv]{spbp}} fit.
#' @param ... arguments inherent from \code{\link[stats]{model.matrix}}.
#' @return The model matrix for the fit.
#' @seealso \code{\link[spsurv]{spbp}}, \code{\link[stats]{model.matrix}}
#' @examples
#' fit <- spbp(Surv(time, status) ~ age + sex, approach = "bayes", data = lung)
#' model.matrix(fit)


model.matrix.spbp <-
  function(spbp, data = eval(spbp$call$data)){
    model.matrix(as.formula(spbp$call$formula), data = data)[, -1]
}
