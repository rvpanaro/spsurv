# Confidence intervals for the regression coefficients

Confidence intervals for the regression coefficients

## Usage

``` r
# S3 method for class 'spbp'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  a fitted model object.

- parm:

  a specification of which parameters are to be given confidence
  intervals: regression coefficient names and/or Bernstein baseline
  names (e.g. `"gamma[1]"`). If missing, all regression coefficients are
  used.

- level:

  the confidence level required.

- ...:

  further arguments passed to parent method

## Value

100(1-alpha) confidence intervals for the requested parameters.
