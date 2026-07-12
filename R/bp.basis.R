#' Bernstein basis polynomials calculations
#' @export
#' @param time a vector of times.
#' @param degree Bernstein polynomial degree
#' @param tau must be greater than times maximum value observed.
#' @return A list containing matrices g and G corresponding BP basis and corresponding tau value used to compute them.
#' @importFrom viridis viridis
#'
bp.basis <- function(time, degree, tau = max(time)) {
  n <- length(time)

  if (sum(time >= 0) != n) stop("time must be a positive vector.")
  if (degree < 0) stop("polynomial degree must be positive.")
  if (!(degree %% 1 == 0)) stop("polynomial degree must be integer.")
  if (tau < max(time)) stop("tau must be greater than the last time.")

  k <- seq_len(degree)
  y <- time / tau

  g <- sapply(k, function(k) {
    dbeta(y, k, degree - k + 1) / tau
  })
  G <- sapply(k, function(k) {
    pbeta(y, k, degree - k + 1)
  })
  list(g = g, G = G, degree = degree, tau = tau)
}

#' @export
#' @method plot spbp
#' @title BP based models plot.
#' @description Plot for a fitted \code{\link[spsurv]{spbp}} model.
#' @param x an object of class `spbp` result of a \code{\link[spsurv]{spbp}} fit.
#' @param main graph title
#' @param graph type of polynomial graph, default is "basis"
#' @param cumulative TRUE for odds and cumulative hazard
#' @param frame graphical parameter; default is FALSE
#' @param lwd graphical parameter; default is 3
#' @param ... further arguments passed to or from other methods
#' @seealso \code{\link[spsurv]{spbp}}.
#' @importFrom graphics curve
#' @examples
#'
#' library("spsurv")
#' data("veteran", package = "survival")
#'
#' fit <- bpph(Surv(time, status) ~ karno + factor(celltype),
#'   data = veteran
#' )
#' plot(fit)
plot.spbp <- function(x, main,
                      graph = c("baseline", "basis"),
                      cumulative = FALSE,
                      frame = FALSE,
                      lwd = 3,
                      ...) {
  graph <- match.arg(graph)
  m <- length(x$bp.param)
  tau_a <- x$tau_a
  tau_b <- x$tau_b
  gamma <- x$bp.param
  pal <- viridis::viridis(m + 1)

  if (missing(main)) {
    if (graph == "basis") {
      main <- if (cumulative) "Cumulative Basis" else "Basis"
    } else if (cumulative) {
      main <- if (x$call$model %in% c("ph", "aft")) {
        "Baseline cumulative hazard"
      } else {
        "Baseline odds"
      }
    } else {
      main <- if (x$call$model %in% c("ph", "aft")) {
        "Baseline hazard"
      } else {
        "Baseline derivative odds"
      }
    }
  }

  portion <- function(x, k) {
    u <- (x - tau_a) / (tau_b - tau_a)
    if (graph == "baseline") {
      if (cumulative) {
        gamma[k] * pbeta(u, k, m - k + 1)
      } else {
        gamma[k] * dbeta(u, k, m - k + 1) / (tau_b - tau_a)
      }
    } else if (cumulative) {
      pbeta(u, k, m - k + 1)
    } else {
      dbeta(u, k, m - k + 1) / (tau_b - tau_a)
    }
  }

  summation <- function(x) {
    vapply(x, function(i) sum(portion(i, seq_len(m))), numeric(1))
  }

  y_max <- if (graph == "basis" && cumulative) {
    1
  } else if (graph == "basis") {
    xs <- seq(tau_a, tau_b, length.out = 500)
    component_max <- max(vapply(seq_len(m), function(k) max(portion(xs, k)), numeric(1)))
    max(component_max, max(summation(xs)))
  } else if (cumulative) {
    max(x$bp.param)
  } else {
    max(x$bp.param) * m / (tau_b - tau_a)
  }

  curve(
    portion(x, 1),
    xlim = c(tau_a, tau_b),
    ylim = c(0, y_max),
    frame = frame,
    lwd = lwd,
    col = pal[1],
    ylab = "BP",
    xlab = "Time",
    main = main,
    ...
  )

  for (k in seq_len(m)) {
    curve(portion(x, k), add = TRUE, col = pal[k + 1], lwd = lwd)
  }

  curve(summation, add = TRUE, lwd = lwd, lty = 2)
}
