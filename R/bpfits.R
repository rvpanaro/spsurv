#' @keywords internal
#' @noRd
.spbp_fit_wrapper <- function(formula, degree = NULL, data, model, approach = c("mle", "bayes"),
                              dist = NULL, baseline = NULL, ..., .mc = NULL) {
  args <- list(
    formula = formula,
    degree = degree,
    data = data,
    model = model,
    approach = match.arg(approach)
  )
  if (!is.null(dist)) {
    args$dist <- dist
  }
  if (!is.null(baseline)) {
    args$baseline <- baseline
  }
  fit <- do.call(spbp.default, c(args, list(...)))
  resolved_call <- fit$call
  mc <- .mc
  if (is.null(mc)) {
    mc <- match.call(expand.dots = TRUE)
  }
  if (is.null(mc$degree)) {
    mc$degree <- NULL
  }
  if (is.null(dist)) {
    mc$dist <- NULL
  }
  if (is.null(baseline)) {
    mc$baseline <- NULL
  }
  fit$call <- mc
  fit$call$approach <- resolved_call$approach
  fit$call$model <- resolved_call$model
  fit
}

#' Bernstein PH Model
#'
#' @export
#' @description Fits the BPPH model to time-to-event data.
#' @param formula a Surv object with time to event observations, right censoring status and explanatory terms.
#' @param degree Bernstein polynomial degree. If omitted, defaults to
#'   \code{ceiling(sqrt(n))} for \code{n = nrow(data)} (see \code{\link{spbp}}).
#' @param data a data.frame object.
#' @param approach Bayesian or maximum likelihood estimation methods, default is approach = "mle".
#' @param dist optional baseline specification; use \code{\link{bernstein}(m)} for the Bernstein polynomial degree.
#' @param baseline optional alias for \code{dist}.
#' @param ... further arguments passed to or from other methods
#' @return An object of class \code{spbp}, including component \code{degree}.
#' @seealso \code{\link[spsurv]{spbp}}, \code{\link[spsurv]{bppo}} and \code{\link[spsurv]{bpaft}} for other BP based models.
#' @examples
#'
#' library("spsurv")
#' data("veteran", package = "survival")
#'
#' fit <- bpph(Surv(time, status) ~ karno + factor(celltype),
#'   data = veteran
#' )
#'
#' summary(fit)
#' @importFrom rstan stan sampling optimizing
#' @importFrom survival Surv frailty


bpph <- function(formula, degree = NULL, data, approach = c("mle", "bayes"),
                 dist = NULL, baseline = NULL, ...) {
  .spbp_fit_wrapper(
    formula = formula,
    degree = degree,
    data = data,
    model = "ph",
    approach = approach,
    dist = dist,
    baseline = baseline,
    ...,
    .mc = match.call(expand.dots = TRUE)
  )
}

#' Bernstein PO Model
#'
#' @export
#' @description Fits the BPPO model to time-to-event data.
#' @param formula a Surv object with time-to-event observations, right censoring status and explanatory terms.
#' @param degree Bernstein polynomial degree. If omitted, defaults to
#'   \code{ceiling(sqrt(n))} for \code{n = nrow(data)} (see \code{\link{spbp}}).
#' @param data a data.frame object.
#' @param approach Bayesian or maximum likelihood estimation methods, default is approach = "mle".
#' @param dist optional baseline specification; use \code{\link{bernstein}(m)} for the Bernstein polynomial degree.
#' @param baseline optional alias for \code{dist}.
#' @param ... further arguments passed to or from other methods
#' @return An object of class \code{spbp}, including component \code{degree}.
#' @seealso \code{\link[spsurv]{spbp}}, \code{\link[spsurv]{bpph}} and \code{\link[spsurv]{bpaft}} for other BP based models.
#' @examples
#'
#' library("spsurv")
#' data("veteran", package = "survival")
#'
#' fit <- bppo(Surv(time, status) ~ karno + celltype,
#'   data = veteran
#' )
#'
#' summary(fit)
#' @importFrom rstan stan sampling optimizing
#' @importFrom survival Surv frailty

bppo <- function(formula, degree = NULL, data, approach = c("mle", "bayes"),
                 dist = NULL, baseline = NULL, ...) {
  .spbp_fit_wrapper(
    formula = formula,
    degree = degree,
    data = data,
    model = "po",
    approach = approach,
    dist = dist,
    baseline = baseline,
    ...,
    .mc = match.call(expand.dots = TRUE)
  )
}

#' Bernstein AFT Model
#'
#' @export
#' @description Fits the BPAFT model to time-to-event data.
#' @param formula a Surv object with time to event observations, right censoring status and explanatory terms.
#' @param degree Bernstein polynomial degree. If omitted, defaults to
#'   \code{ceiling(sqrt(n))} for \code{n = nrow(data)} (see \code{\link{spbp}}).
#' @param data a data.frame object.
#' @param approach Bayesian or maximum likelihood estimation methods, default is approach = "mle".
#' @param dist optional baseline specification; use \code{\link{bernstein}(m)} for the Bernstein polynomial degree.
#' @param baseline optional alias for \code{dist}.
#' @param ... further arguments passed to or from other methods
#' @return An object of class \code{spbp}, including component \code{degree}.
#' @seealso \code{\link[spsurv]{spbp}}, \code{\link[spsurv]{bpph}} and \code{\link[spsurv]{bppo}} for other BP based models.
#' @examples
#'
#' library("spsurv")
#' data("veteran", package = "survival")
#'
#' fit <- bpaft(Surv(time, status) ~ karno + celltype,
#'   data = veteran
#' )
#'
#' summary(fit)
#' @importFrom rstan stan sampling optimizing
#' @importFrom survival Surv frailty

bpaft <- function(formula, degree = NULL, data, approach = c("mle", "bayes"),
                  dist = NULL, baseline = NULL, ...) {
  .spbp_fit_wrapper(
    formula = formula,
    degree = degree,
    data = data,
    model = "aft",
    approach = approach,
    dist = dist,
    baseline = baseline,
    ...,
    .mc = match.call(expand.dots = TRUE)
  )
}
