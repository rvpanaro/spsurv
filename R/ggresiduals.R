#' ggplot2 residual diagnostic plots for spbp models
#'
#' Convenience wrapper around \code{\link{residuals.spbp}} that returns a
#' \pkg{ggplot2} object for martingale, deviance, or Cox-Snell residuals.
#'
#' @param object A fitted \code{"spbp"} object.
#' @param type Residual type: \code{"martingale"} (default), \code{"deviance"},
#'   \code{"cox-snell"}, or \code{"coxsnell"}.
#' @param against What to plot on the x-axis: \code{"fitted"} (default, Cox-Snell
#'   cumulative hazard) or \code{"index"} (observation index).
#' @param ... Further arguments passed to \code{\link{residuals.spbp}}.
#' @return A \code{ggplot} object (requires \pkg{ggplot2}).
#' @export
#' @examples
#' \dontrun{
#' library(spsurv)
#' library(ggplot2)
#' data(veteran, package = "survival")
#' fit <- bpph(Surv(time, status) ~ karno, data = veteran, degree = 4)
#' ggresiduals(fit, type = "martingale")
#' }
ggresiduals <- function(object,
                        type = c("martingale", "deviance", "cox-snell", "coxsnell"),
                        against = c("fitted", "index"),
                        ...) {
  if (!inherits(object, "spbp")) {
    stop("Object must inherit from class 'spbp'.", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for ggresiduals().", call. = FALSE)
  }

  type <- match.arg(type)
  against <- match.arg(against)
  if (identical(type, "coxsnell")) {
    type <- "cox-snell"
  }

  res <- residuals(object, type = type, ...)
  if (against == "fitted") {
    xvar <- residuals(object, type = "cox-snell", ...)
    xlab <- "Fitted cumulative hazard (Cox-Snell)"
  } else {
    xvar <- seq_along(res)
    xlab <- "Observation index"
  }

  df <- data.frame(x = xvar, y = res)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.8) +
    ggplot2::labs(
      x = xlab,
      y = paste(type, "residual"),
      title = paste(type, "residuals")
    ) +
    ggplot2::theme_bw()

  if (against == "fitted") {
    p <- p +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      ggplot2::geom_smooth(
        mapping = ggplot2::aes(group = 1),
        method = "loess",
        se = FALSE,
        color = "black",
        linewidth = 0.5
      )
  }

  p
}
