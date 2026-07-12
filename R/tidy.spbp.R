#' Tidy summary for fitted spbp models
#'
#' @description
#' Create a broom/generics-style tidy data frame with one row per model term.
#'
#' @param x A fitted \code{"spbp"} object.
#' @param conf.int Logical; include confidence/credible interval columns.
#' @param conf.level Interval coverage level.
#' @param exponentiate Logical; if \code{TRUE}, exponentiate the estimate and
#'   interval columns.
#' @param ... Currently unused.
#'
#' @return A \code{data.frame} with columns \code{term}, \code{estimate},
#'   \code{std.error}, and when available \code{statistic}, \code{p.value},
#'   \code{conf.low}, \code{conf.high}. Columns that are entirely
#'   \code{NA} are omitted.
#' @export
#' @method tidy spbp
#' @importFrom generics tidy
tidy.spbp <- function(x, conf.int = FALSE, conf.level = 0.95,
                      exponentiate = FALSE, ...) {
  if (is.null(x$coefficients)) {
    return(data.frame(
      term = character(0),
      estimate = numeric(0),
      std.error = numeric(0),
      statistic = numeric(0),
      p.value = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  sm <- summary(x, interval = conf.level)
  coef_mat <- sm$coefficients

  if (identical(x$call$approach, "mle")) {
    out <- data.frame(
      term = rownames(coef_mat),
      estimate = coef_mat[, "coef"],
      std.error = coef_mat[, "se(coef)"],
      statistic = coef_mat[, "z"],
      p.value = coef_mat[, "Pr(>|z|)"],
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  } else {
    out <- data.frame(
      term = rownames(coef_mat),
      estimate = coef_mat[, "mean(coef)"],
      std.error = coef_mat[, "sd(coef)"],
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }

  if (isTRUE(conf.int)) {
    if (identical(x$call$approach, "mle")) {
      ci <- stats::confint(x, level = conf.level)
    } else {
      ci <- credint(x, prob = conf.level, type = "HPD")
    }
    ci <- as.matrix(ci)
    ci <- ci[match(out$term, rownames(ci)), , drop = FALSE]
    out$conf.low <- ci[, 1]
    out$conf.high <- ci[, 2]
  }

  if (isTRUE(exponentiate)) {
    out$estimate <- exp(out$estimate)
    if (isTRUE(conf.int)) {
      out$conf.low <- exp(out$conf.low)
      out$conf.high <- exp(out$conf.high)
    }
  }

  .spbp_drop_all_na_cols(out)
}
