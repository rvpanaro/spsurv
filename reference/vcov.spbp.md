# Covariance of the regression coefficients

Uses block-wise inversion of the negative Hessian, with a clear split
between the regression coefficients (beta) and the Bernstein polynomial
coefficients (gamma).

## Usage

``` r
# S3 method for class 'spbp'
vcov(object, bp.param = FALSE, ...)
```

## Arguments

- object:

  an object of the class spbp

- bp.param:

  return Bernstein Polynomial variance.

- ...:

  arguments passed to parent method.

## Value

the variance-covariance matrix associated with the regression
coefficients.
