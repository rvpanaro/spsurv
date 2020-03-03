#' Bernstein polynomial based proportional hazards model
#'
#' @export
#' @description Fits the BPPH model to time-to-event data.
#' @param formula a Surv object with time to event observations, right censoring status and explanatory terms.
#' @param degree Bernstein polynomial degree.
#' @param data a data.frame object.
#' @param approach Bayesian or maximum likelihood estimation methods, default is approach = "mle".
#' @param priors Prior choice for the Bayesian analysis, see \code{\link[spsurv]{spbp}}.
#' @return An object of class `spbp`.
#' @seealso \code{\link[spsurv]{spbp}}, \code{\link[spsurv]{bppo}} and \code{\link[spsurv]{bpaft}} for other BP based models.
#' @examples
#'
#' library("KMsurv")
#' data("larynx")
#'
#' library("spsurv")
#' fit <- bpph(Surv(time, status) ~ age + factor(stage),
#' data = larynx)
#'
#' summary(fit)
#'
  #' @references Osman, M. and Ghosh, S. K. (2012), “Nonparametric regression models for right-censoreddata using Bernstein polynomials,”Computational Statistics & Data Analysis, 56, 559–573.
#' @importFrom rstan stan sampling optimizing
#' @importFrom survival Surv frailty


bpph <- function(formula, degree, data,  approach = c("mle", "bayes"), ...){
  fit <- spbp.default(formula = formula,
               degree = degree,
               data = data,
               model = "ph",
               approach = match.arg(approach),
               ...)
  fit$call$formula <- match.call()$formula
  fit$call$data <- match.call()$data
  return(fit)
}

#' Bernstein polynomial based proportional odds model
#'
#' @export
#' @description Fits the BPPO model to time-to-event data.
#' @param formula a Surv object with time-to-event observations, right censoring status and explanatory terms.
#' @param degree Bernstein polynomial degree.
#' @param data a data.frame object.
#' @param approach Bayesian or maximum likelihood estimation methods, default is approach = "mle".
#' @param priors Prior choice for the Bayesian analysis, see \code{\link[spsurv]{spbp}}.
#' @return An object of class `spbp`.
#' @seealso \code{\link[spsurv]{spbp}}, \code{\link[spsurv]{bpph}} and \code{\link[spsurv]{bpaft}} for other BP based models.
#' @examples
#'
#' library("survival")
#' data("veteran")
#'
#' library("spsurv")
#' fit <- bppo(Surv(time, status) ~ karno + celltype,
#' data = veteran)
#'
#' summary(fit)
#'
#' @references Osman, M. and Ghosh, S. K. (2012), “Nonparametric regression models for right-censoreddata using Bernstein polynomials,”Computational Statistics & Data Analysis, 56, 559–573.
#' @importFrom rstan stan sampling optimizing
#' @importFrom survival Surv frailty

bppo <- function(formula, degree, data,  approach = c("mle", "bayes"), ...){
  fit <- spbp.default(formula = formula,
               degree = degree,
               data = data,
               model = "po",
               approach = match.arg(approach),
               ...)
  fit$call$formula <- match.call()$formula
  fit$call$data <- match.call()$data
  return(fit)
}

#' Bernstein polynomial based Accelerated Failure Time
#'
#' @export
#' @description Fits the BPAFT model to time-to-event data.
#' @param formula a Surv object with time to event observations, right censoring status and explanatory terms.
#' @param degree Bernstein polynomial degree.
#' @param data a data.frame object.
#' @param approach Bayesian or maximum likelihood estimation methods, default is approach = "mle".
#' @param priors Prior choice for the Bayesian analysis, see \code{\link[spsurv]{spbp}}.
#' @return An object of class `spbp`.
#' @seealso \code{\link[spsurv]{spbp}}, \code{\link[spsurv]{bpph}} and \code{\link[spsurv]{bppo}} for other BP based models.
#' @examples
#'
#' library("survival")
#' data("veteran")
#'
#' library("spsurv")
#' fit <- bpaft(Surv(time, status) ~ karno + celltype,
#' data = veteran)
#'
#' summary(fit)
#'
#' @references Osman, M. and Ghosh, S. K. (2012), “Nonparametric regression models for right-censoreddata using Bernstein polynomials,”Computational Statistics & Data Analysis, 56, 559–573.
#' @importFrom rstan stan sampling optimizing
#' @importFrom survival Surv frailty

bpaft <- function(formula, degree, data, approach = c("mle", "bayes"), ...){
  fit <- spbp.default(formula =  formula,
               degree = degree,
               data = data,
               model = "aft",
               approach = match.arg(approach),
               ...)
  fit$call$formula <- match.call()$formula
  fit$call$data <- match.call()$data
  return(fit)
}
