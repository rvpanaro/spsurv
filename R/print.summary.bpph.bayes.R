#' Bernstein Polynomial Based Regression Object Summary BPPH Bayes
#'
#' @export
#' @param ... further arguments passed to or from other methods
#' @method print summary.bpph.bayes
#' @return none

print.summary.bpph.bayes <-
  function(...) {
    print.summary.spbp.bayes(...)
  }
