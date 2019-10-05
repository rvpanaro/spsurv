#' Bernstein Polynomial based Regression Object
#'
#' @export
#' @method print summary.bpaft.bayes
#'

print.summary.bpaft.bayes <-
  function(...){
  cat("Bayesian Bernstein Polynomial based Accelerated Failure Time model\n")
  print.summary.spbp.bayes(...)
}
