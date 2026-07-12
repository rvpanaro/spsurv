#' @keywords internal
#' @noRd
.spsurv_parsnip_warn_weights <- function(case_weights) {
  if (!is.null(case_weights)) {
    warning("Case weights are not supported for spsurv parsnip engines.", call. = FALSE)
  }
}

#' @keywords internal
#' @noRd
.spsurv_parsnip_fit <- function(formula, data, model = c("ph", "po", "aft"),
                                case_weights = NULL, degree = NULL,
                                approach = c("mle", "bayes"),
                                scale = FALSE, init = 0, verbose = FALSE,
                                ...) {
  .spsurv_parsnip_warn_weights(case_weights)
  model <- match.arg(model)
  approach <- match.arg(approach)
  fit_fun <- switch(
    model,
    ph = bpph,
    po = bppo,
    aft = bpaft
  )
  fit_fun(
    formula = formula,
    data = data,
    degree = degree,
    approach = approach,
    scale = scale,
    init = init,
    verbose = verbose,
    ...
  )
}

#' Fit proportional hazards model (spsurv parsnip engine)
#' @keywords internal
#' @export
spsurv_fit_proportional_hazards <- function(formula, data, case_weights = NULL, ...) {
  .spsurv_parsnip_fit(formula, data, model = "ph", case_weights = case_weights, ...)
}

#' Fit proportional odds model (spsurv parsnip engine)
#' @keywords internal
#' @export
spsurv_fit_proportional_odds <- function(formula, data, case_weights = NULL, ...) {
  .spsurv_parsnip_fit(formula, data, model = "po", case_weights = case_weights, ...)
}

#' Fit AFT model (spsurv parsnip engine)
#' @keywords internal
#' @export
spsurv_fit_survival_reg <- function(formula, data, case_weights = NULL, ...) {
  .spsurv_parsnip_fit(formula, data, model = "aft", case_weights = case_weights, ...)
}
