#' Bernstein Polynomial based Regression Object
#'
#' @export
#' @param ... further arguments passed to or from other methods
#' @method print summary.bppo.mle
#'
print.summary.bppo.mle <-
  function(...){
  cat("Bernstein Polynomial based Proportional Odds model\n")
  print.summary.spbp.mle(...)
}
