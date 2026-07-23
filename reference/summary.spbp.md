# Bernstein Polynomial Based Regression Object Summary

Bernstein Polynomial Based Regression Object Summary

## Usage

``` r
# S3 method for class 'spbp'
summary(
  object,
  interval = 0.95,
  compact = TRUE,
  show_call = TRUE,
  show_intervals = TRUE,
  show_baseline = FALSE,
  mle_test = c("lr", "wald"),
  bayes_criterion = c("waic", "dic", "lpml"),
  ...
)
```

## Arguments

- object:

  an object of class spbp

- interval:

  interval coverage (confidence or credibility)

- compact:

  logical; if TRUE, print.summary methods show essential output only.

- show_call:

  logical; include the model call in printed summary.

- show_intervals:

  logical; include interval table in printed summary.

- show_baseline:

  logical; include Bernstein baseline (`gamma`) Wald table in printed
  summary (MLE only).

- mle_test:

  global test to print for MLE summaries: "lr" or "wald".

- bayes_criterion:

  global criterion to print for Bayesian summaries: "waic", "dic", or
  "lpml".

- ...:

  further arguments passed to or from other methods

## Value

An object of class analogous to for e.g. 'summary.bppo.bayes'.
