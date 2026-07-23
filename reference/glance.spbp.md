# Glance at a fitted spbp model

Return a one-row broom/generics-style summary of model-level statistics.

## Usage

``` r
# S3 method for class 'spbp'
glance(x, ...)
```

## Arguments

- x:

  A fitted `"spbp"` object.

- ...:

  Currently unused.

## Value

A one-row `data.frame` with model-level statistics. Columns that are
entirely `NA` for the fitted object are omitted.
