# Confidence intervals for the regression coefficients

Confidence intervals for the regression coefficients

## Usage

``` r
# S3 method for class 'spbp'
confint(object, parm = names(coef(object)), level = 0.95, ...)
```

## Arguments

- object:

  a fitted model object.

- parm:

  a specification of which parameters are to be given confidence
  intervals, either a vector of numbers or a vector of names. If
  missing, all parameters are considered.

- level:

  the confidence level required.

- ...:

  further arguments passed to parent method

## Value

100(1-alpha) confidence intervals for the regression coefficients
