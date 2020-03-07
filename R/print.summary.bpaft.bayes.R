#' Bernstein Polynomial based Regression Object
#'
#' @export
#' @param ... further arguments passed to or from other methods
#' @method print summary.bpaft.bayes
#'

print.summary.bpaft.bayes <-
  function(...){
  cat("Bayesian Bernstein Polynomial based Accelerated Failure Time model\n")
  print.summary.spbp.bayes(...)
}
