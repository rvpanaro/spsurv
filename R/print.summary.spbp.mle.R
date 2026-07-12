#' Bernstein Polynomial Based Regression Object Summary MLE
#'
#' @export
#' @param x a summary.spbp.mle object
#' @param digits number of digits to display.
#' @param signif.stars see \code{\link{getOption}}
#' @param ... further arguments passed to or from other methods
#' @method print summary.spbp.mle
#' @return none

print.summary.spbp.mle <-
  function(x, digits = 2,
           signif.stars = getOption("show.signif.stars"), ...) {
    .print.summary.spbp(
      x,
      digits = digits,
      signif.stars = signif.stars,
      approach = "mle",
      ...
    )
  }
