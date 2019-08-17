traceplot.spbp <- function(spbp, pars = 'beta', ...){
  if(spbp$approach != 1)
    stop("Invalid input, approach should be 'bayes'.")
  rstan::traceplot(spbp$stan, pars = pars, ...)
}
#
# traceplot.bppo.bayes <- function( ...){
#   traceplot.spbp(...)
# }
#
# traceplot.bpph.bayes <- function(...){
#   traceplot.spbp(...)
# }
#
#

survfit.spbp <- function(spbp){
  print('Echo')
}
