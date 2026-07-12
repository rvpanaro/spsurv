#' Augment data with spbp model information
#'
#' @description
#' Add fitted values, residuals, and optional survival predictions to the
#' training (or supplied) data.
#'
#' @param x A fitted \code{"spbp"} object.
#' @param data Optional data frame; defaults to the training data when available.
#' @param eval_time Optional numeric vector of times for nested survival
#'   predictions (same structure as \code{predict(x, type = "survival")}).
#' @param type Residual type passed to \code{\link{residuals.spbp}}.
#' @param ... Not used.
#' @return A \code{data.frame} with original rows plus \code{.residual} and,
#'   when \code{eval_time} is set, a list-column \code{.pred}.
#' @export
#' @method augment spbp
#' @importFrom generics augment
augment.spbp <- function(x, data = NULL, eval_time = NULL,
                         type = c("martingale", "deviance", "cox-snell", "coxsnell"), ...) {
  if (is.null(data)) {
    data <- x$data
    if (is.null(data)) {
      data <- .spbp_training_data(x)
    }
  }
  out <- as.data.frame(data)
  res <- residuals(x, type = type)
  if (length(res) == nrow(out)) {
    out$.residual <- res
  } else if (!is.null(names(res)) && !is.null(rownames(out))) {
    matched <- res[rownames(out)]
    if (sum(!is.na(matched)) > 0L) {
      out$.residual <- matched
    }
  }
  if (!is.null(eval_time)) {
    pred <- .spbp_predict_survival_censored(x, newdata = data, eval_time = eval_time)
    out$.pred <- pred$.pred
  }
  out
}
