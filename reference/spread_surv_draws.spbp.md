# Posterior survival curves in long tidy format

Posterior survival curves in long tidy format

## Usage

``` r
spread_surv_draws.spbp(object, times, newdata = NULL, ...)
```

## Arguments

- object:

  A Bayesian `"spbp"` fit.

- times:

  Numeric vector of evaluation times.

- newdata:

  Optional covariate profiles.

- ...:

  Not used.

## Value

A `data.frame` with `.chain`, `.iteration`, `.draw`, `time`, `surv`, and
optional covariate columns.
