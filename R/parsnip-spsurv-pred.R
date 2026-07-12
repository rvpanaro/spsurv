#' @keywords internal
#' @noRd
.spsurv_parsnip_extract_fit <- function(object) {
  if (inherits(object, "spbp")) {
    return(object)
  }
  if (is.list(object) && inherits(object$fit, "spbp")) {
    return(object$fit)
  }
  stop("Expected a parsnip fit containing an 'spbp' object.", call. = FALSE)
}

#' Predict survival probabilities (spsurv parsnip engine)
#' @keywords internal
#' @export
spsurv_pred_survival <- function(object, new_data, eval_time, ...) {
  fit <- .spsurv_parsnip_extract_fit(object)
  predict(fit, newdata = new_data, type = "survival", eval_time = eval_time)
}

#' Predict median event time (spsurv parsnip engine)
#' @keywords internal
#' @export
spsurv_pred_time <- function(object, new_data, ...) {
  fit <- .spsurv_parsnip_extract_fit(object)
  predict(fit, newdata = new_data, type = "time")
}

#' Predict linear predictor (spsurv parsnip engine)
#' @keywords internal
#' @export
spsurv_pred_linear_pred <- function(object, new_data, ...) {
  fit <- .spsurv_parsnip_extract_fit(object)
  predict(fit, newdata = new_data, type = "linear_pred")
}
