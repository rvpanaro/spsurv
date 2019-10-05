#' Bernstein Polynomial based Regression Object
#'
#' @export
#' @method print summary.bpaft.mle
#'

print.summary.bpaft.mle <-
  function(...){
  cat("Bernstein Polynomial based Accelerated Failure Time model\n")
  print.summary.spbp.mle(...)
}
