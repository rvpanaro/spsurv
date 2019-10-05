#' Bernstein Polynomial based Regression Object
#'
#' @export
#' @method print summary.bpph.mle
#'

print.summary.bpph.mle <-
  function(...){
  cat("Bernstein Polynomial based Proportional Hazards model\n")
  print.summary.spbp.mle(...)
}
