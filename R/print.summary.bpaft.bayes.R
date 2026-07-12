#' Bernstein Polynomial Based Regression Object Summary BPAFT Bayes
#'
#' @export
#' @param ... further arguments passed to or from other methods
#' @method print summary.bpaft.bayes
#' @return none

print.summary.bpaft.bayes <-
  function(...) {
    print.summary.spbp.bayes(...)
  }
