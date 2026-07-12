# Quick Bayesian fits for tests (short Stan runs; warnings are non-deterministic).
quick_bayes <- function(fit_fun, ...) {
  suppressWarnings(fit_fun(...))
}
