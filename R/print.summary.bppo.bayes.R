#' Bernstein Polynomial based Regression Object
#'
#' @export
#' @param ... further arguments passed to or from other methods
#' @method print summary.bppo.bayes
#'

print.summary.bppo.bayes <-
  function(...){
  cat("Bayesian Bernstein Polynomial based Proportional Odds model\n")
  print.summary.spbp.bayes(...)
}
