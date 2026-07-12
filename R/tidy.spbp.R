#' Tidy summary for fitted spbp models
#'
#' @description
#' Create a broom/generics-style tidy data frame with one row per model term.
#'
#' @param x A fitted \code{"spbp"} object.
#' @param conf.int Logical; include confidence/credible interval columns.
#' @param conf.level Interval coverage level.
#' @param exponentiate Logical; if \code{TRUE}, exponentiate the estimate and
#'   interval columns (regression coefficients only).
#' @param component Which parameters to return: \code{"coef"} (default),
#'   \code{"baseline"} (Bernstein \code{gamma} weights), or \code{"all"}.
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
                      exponentiate = FALSE,
                      component = c("coef", "baseline", "all"), ...) {
  component <- match.arg(component)

  if (component %in% c("baseline", "all") && identical(x$call$approach, "mle")) {
    base_tab <- .spbp_baseline_mle_table(x)
    base_out <- data.frame(
      term = rownames(base_tab),
      component = "baseline",
      estimate = base_tab[, "gamma"],
      std.error = base_tab[, "se_log_gamma"],
      statistic = base_tab[, "z"],
      p.value = base_tab[, "Pr(>|z|)"],
      row.names = NULL,
      stringsAsFactors = FALSE
    )
    if (isTRUE(conf.int)) {
      ci <- confint(x, parm = rownames(base_tab), level = conf.level)
      ci <- as.matrix(ci)
      base_out$conf.low <- ci[, 1]
      base_out$conf.high <- ci[, 2]
    }
    base_out <- .spbp_drop_all_na_cols(base_out)
  } else if (component %in% c("baseline", "all")) {
    base_out <- NULL
  } else {
    base_out <- NULL
  }

  if (component %in% c("coef", "all")) {
    if (is.null(x$coefficients)) {
      coef_out <- data.frame(
        term = character(0),
        component = character(0),
        estimate = numeric(0),
        std.error = numeric(0),
        statistic = numeric(0),
        p.value = numeric(0),
        stringsAsFactors = FALSE
      )
    } else {
      sm <- summary(x, interval = conf.level)
      coef_mat <- sm$coefficients

      if (identical(x$call$approach, "mle")) {
        coef_out <- data.frame(
          term = rownames(coef_mat),
          component = "coef",
          estimate = coef_mat[, "coef"],
          std.error = coef_mat[, "se(coef)"],
          statistic = coef_mat[, "z"],
          p.value = coef_mat[, "Pr(>|z|)"],
          row.names = NULL,
          stringsAsFactors = FALSE
        )
      } else {
        coef_out <- data.frame(
          term = rownames(coef_mat),
          component = "coef",
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
        ci <- ci[match(coef_out$term, rownames(ci)), , drop = FALSE]
        coef_out$conf.low <- ci[, 1]
        coef_out$conf.high <- ci[, 2]
      }

      if (isTRUE(exponentiate)) {
        coef_out$estimate <- exp(coef_out$estimate)
        if (isTRUE(conf.int)) {
          coef_out$conf.low <- exp(coef_out$conf.low)
          coef_out$conf.high <- exp(coef_out$conf.high)
        }
      }

      coef_out <- .spbp_drop_all_na_cols(coef_out)
    }
  } else {
    coef_out <- NULL
  }

  if (component == "coef") {
    return(coef_out)
  }
  if (component == "baseline") {
    if (is.null(base_out)) {
      stop(
        "Baseline tidy output is only available for MLE fits.",
        call. = FALSE
      )
    }
    return(base_out)
  }

  out <- rbind(coef_out, base_out)
  row.names(out) <- NULL
  out
}
