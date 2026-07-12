#' Bernstein Polynomial Based Regression Object Summary Bayes
#'
#' @export
#' @param x a summary.spbp.bayes object
#' @param digits number of digits to display.
#' @param signif.stars see \code{\link{getOption}}
#' @param ... further arguments passed to or from other methods
#' @method print summary.spbp.bayes
#' @return none

print.summary.spbp.bayes <-
  function(x, digits = 2,
           signif.stars = getOption("show.signif.stars"), ...) {
    .print.summary.spbp(
      x,
      digits = digits,
      signif.stars = signif.stars,
      approach = "bayes",
      ...
    )
  }
