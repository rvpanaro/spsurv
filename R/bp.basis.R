#' Bernstein basis polynomials calculations
#' @export
#' @param time a vector of times.
#' @param degree Bernstein polynomial degree
#' @param tau must be greater than times maximum value observed.
#' @return A list containing matrices g and G corresponding BP basis and corresponding tau value used to compute them.

bp.basis <- function(time, degree, tau = max(time)) {
  n <- length(time)

  ## error handling
  if (sum(time >= 0) != n) {
    stop("time must be a positive vector.")
  }

  if (degree < 0) {
    stop("polynomial degree must be positive.")
  }

  if (!(degree %% 1 == 0)) {
    stop("polynomial degree must be integer.")
  }

  if (tau < max(time)) {
    stop("tau must be greater than the last time.")
  }

  k <- 1:degree
  g <- matrix(NA, n, degree)
  G <- matrix(NA, n, degree)
  y <- time / tau

  g <- sapply(k, function(k) {
    dbeta(y, k, degree - k + 1) / tau
  })
  G <- sapply(k, function(k) pbeta(y, k, degree - k + 1))

  # Equivalent to
  # for (i in 1:n){
  #   for(k in 1:degree){
  #     g[i,k] <- dbeta(y[i], k, degree - k + 1) / tau
  #     G[i,k] <- pbeta(y[i], k, degree - k + 1)
  #   }
  # }
  return(list(g = g, G = G, degree = degree, tau = tau))
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
#' @examples
#'
#' library("spsurv")
#' data("veteran")
#'
#' fit <- bpph(Surv(time, status) ~ karno + factor(celltype),
#'   data = veteran
#' )
#' plot(fit)
plot.spbp <- function(x, main,
                      graph = c("baseline", "basis"),
                      cumulative = F,
                      frame = F,
                      lwd = 3,
                      ...) {
  graph <- match.arg(graph)
  m <- length(x$bp.param)
  tau_a <- x$tau_a
  tau_b <- x$tau_b
  approach <- x$call$approach
  gamma <- x$bp.param

  if (missing(main)) {
    if (graph == "basis") {
      if (cumulative) {
        main <- "Cumulative Basis"
      } else {
        main <- "Basis"
      }
    } else if (cumulative) {
      if (x$call$model %in% c("ph", "aft")) {
        main <- "Baseline cumulative hazard"
      } else {
        main <- "Baseline odds"
      }
    } else {
      if (x$call$model %in% c("ph", "aft")) {
        main <- "Baseline hazard"
      } else {
        main <- "Baseline derivative odds"
      }
    }
  }

  portion <- function(x, k) {
    if (graph == "baseline") {
      if (cumulative) {
        gamma[k] * pbeta((x - tau_a) / (tau_b - tau_a), k, m - k + 1)
      } else {
        gamma[k] * dbeta((x - tau_a) / (tau_b - tau_a), k, m - k + 1) / (tau_b - tau_a)
      }
    } else {
      if (cumulative) {
        pbeta((x - tau_a) / (tau_b - tau_a), k, m - k + 1)
      } else {
        dbeta((x - tau_a) / (tau_b - tau_a), k, m - k + 1) / (tau_b - tau_a)
      }
    }
  }

  summation <- function(x) {
    sapply(x, function(i) {
      sum(portion(i, 1:m))
    })
  }

  curve(portion(x, 1),
    xlim = c(tau_a, tau_b),
    ylim = c(0, cumulative * ifelse(graph == "basis", 1, max(x$bp.param)) + (!cumulative) * ifelse(graph == "basis", 1, max(x$bp.param)) * m / (tau_b - tau_a)),
    frame = frame, lwd = lwd,
    col = viridis::viridis(m + 1)[1], ylab = "BP", xlab = "Time",
    main = main, ...
  )

  for (k in 1:m) {
    curve(portion(x, k), add = T, col = viridis::viridis(m + 1)[k + 1], lwd = lwd)
  }

  curve(summation, add = T, lwd = lwd, lty = 2)
}
