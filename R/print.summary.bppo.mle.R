#' Bernstein Polynomial based Regression Object
#'
#' @export
#' @method print summary.bppo.mle
#'
print.summary.bppo.mle <-
  function(...){
  cat("Bernstein Polynomial based Proportional Odds model\n")
  print.summary.spbp.mle(...)
}
