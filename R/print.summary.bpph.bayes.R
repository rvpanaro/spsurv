#' Bernstein Polynomial based Regression Object
#'
#' @export
#' @method print summary.bpph.bayes
#'
print.summary.bpph.bayes <-
  function(...){
  cat("Bayesian Bernstein Polynomial based Proportional Hazards model\n")
  print.summary.spbp.bayes(...)
}
