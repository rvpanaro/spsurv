# Proportional odds regression model specification

Proportional odds survival regression with a Bernstein polynomial
baseline
([`bppo`](https://rvpanaro.github.io/spsurv/reference/bppo.md)). This
model type is registered by spsurv for use with tidymodels; it is not
part of the censored package.

## Usage

``` r
proportional_odds(mode = "censored regression", engine = "spsurv")
```

## Arguments

- mode:

  A single character string for the prediction outcome mode. The only
  possible value for this model is "censored regression".

- engine:

  A single character string specifying what computational engine to use
  for fitting.

## Value

A model specification.
