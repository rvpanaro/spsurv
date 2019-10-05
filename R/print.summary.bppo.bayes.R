#' Bernstein Polynomial based Regression Object
#'
#' @export
#' @method print summary.bppo.bayes
#'

print.summary.bppo.bayes <-
  function(...){
  cat("Bayesian Bernstein Polynomial based Proportional Odds model\n")
  print.summary.spbp.bayes(...)
}
