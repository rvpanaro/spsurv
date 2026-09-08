# Rank fitted spbp models by AIC

Convenience wrapper around [`AIC`](https://rdrr.io/r/stats/AIC.html)
that returns models sorted from lowest to highest AIC.

## Usage

``` r
rank_models(...)
```

## Arguments

- ...:

  One or more fitted `"spbp"` objects (MLE).

## Value

A `data.frame` with columns `fit`, `model`, `degree`, `aic`, and
`npars`.

## Examples

``` r
if (FALSE) { # \dontrun{
library(spsurv)
data(veteran, package = "survival")
f1 <- bpph(Surv(time, status) ~ 1, data = veteran, degree = 4)
f2 <- bpph(Surv(time, status) ~ karno, data = veteran, degree = 4)
rank_models(f1, f2)
} # }
```
