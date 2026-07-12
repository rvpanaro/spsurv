#' Log-likelihood for fitted spbp models
#' @param object A fitted \code{"spbp"} object.
#' @param ... Not used.
#' @export
#' @method logLik spbp
logLik.spbp <- function(object, ...) {
  if (!inherits(object, "spbp")) {
    stop("Object must inherit from class 'spbp'.", call. = FALSE)
  }
  if (!identical(object$call$approach, "mle")) {
    warning(
      "logLik is only defined for MLE fits (approach = 'mle').",
      call. = FALSE
    )
    return(NULL)
  }
  ll <- .spbp_model_loglik(object)
  if (!is.finite(ll)) {
    return(NULL)
  }
  structure(
    ll,
    class = "logLik",
    df = .spbp_nparams(object),
    nobs = object$n
  )
}

#' Akaike information criterion for fitted spbp models
#' @param object A fitted \code{"spbp"} object.
#' @param ... Additional fitted models for comparison.
#' @param k Penalty per parameter (default 2).
#' @export
#' @method AIC spbp
AIC.spbp <- function(object, ..., k = 2) {
  if (!missing(k) && k != 2) {
    warning("Only k = 2 (AIC) is supported for spbp objects.", call. = FALSE)
  }
  dots <- list(...)
  if (length(dots) == 0L) {
    .spbp_check_mle(object)
    return(.spbp_aic(object))
  }

  fits <- c(list(object), dots)
  for (fit in fits) {
    .spbp_check_mle(fit)
  }

  nms <- names(fits)
  if (is.null(nms)) {
    nms <- as.character(seq_along(fits))
  }

  data.frame(
    fit = nms,
    aic = vapply(fits, .spbp_aic, numeric(1)),
    npars = vapply(fits, .spbp_nparams, integer(1)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Extract AIC from a fitted spbp model
#' @param fit A fitted \code{"spbp"} object.
#' @param scale Not used.
#' @param k Penalty per parameter (default 2).
#' @param ... Not used.
#' @export
#' @method extractAIC spbp
extractAIC.spbp <- function(fit, scale = 0, k = 2, ...) {
  .spbp_check_mle(fit)
  c(.spbp_nparams(fit), .spbp_aic(fit))
}

#' Analysis of deviance table for nested spbp models
#'
#' Compare nested MLE fits. With one model, builds a sequential table by refitting
#' intercept-only and nested submodels along the formula term order (same idea as
#' \code{anova(fit)} for \code{\link[survival]{survreg}}). With multiple models,
#' compares each consecutive pair (same convention as
#' \code{anova(fit0, fit1)} for \code{survreg}).
#'
#' @param object A fitted \code{"spbp"} object.
#' @param ... Additional nested fitted models (MLE, same family recommended).
#' @param test Which test to report (only \code{"Chisq"} is supported).
#' @return An \code{"anova"} object (also a \code{data.frame}). Pairwise
#'   comparisons use \code{Terms}, \code{Resid. Df}, \code{-2*LL}, \code{Test},
#'   \code{Df}, \code{Deviance}, and \code{Pr(>Chi)} columns (as in
#'   \code{survival::survreg}). Sequential single-model tables use
#'   \code{Df}, \code{Deviance}, \code{Resid. Df}, \code{-2*LL}, and
#'   \code{Pr(>Chi)} with row names \code{NULL} and term labels.
#' @export
#' @method anova spbp
anova.spbp <- function(object, ..., test = "Chisq") {
  if (!identical(test, "Chisq")) {
    stop("Only test = 'Chisq' is supported.", call. = FALSE)
  }

  dots <- list(...)
  if (length(dots) == 0L) {
    .spbp_check_mle(object)
    return(.spbp_anova_sequential(object))
  }

  fits <- c(list(object), dots)
  for (fit in fits) {
    .spbp_check_mle(fit)
  }

  models <- vapply(fits, function(f) f$call$model, character(1))
  if (length(unique(models)) > 1L) {
    warning(
      "Models have different 'model' families; LR comparisons may not be valid.",
      call. = FALSE
    )
  }

  degrees <- vapply(fits, function(f) f$degree, integer(1))
  if (length(unique(degrees)) > 1L) {
    warning(
      "Models have different Bernstein degrees; LR comparisons may not be valid.",
      call. = FALSE
    )
  }

  n <- length(fits)
  tab <- .spbp_anova_pairwise(fits)
  .spbp_as_anova_table(tab, .spbp_anova_heading(object))
}
