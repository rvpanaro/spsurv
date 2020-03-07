#' Bernstein Polynomial based Regression Object
#'
#' @export
#' @param ... further arguments passed to or from other methods
#' @method print summary.bpaft.mle
#'

print.summary.bpaft.mle <-
  function(...){
  cat("Bernstein Polynomial based Accelerated Failure Time model\n")
  print.summary.spbp.mle(...)
}
