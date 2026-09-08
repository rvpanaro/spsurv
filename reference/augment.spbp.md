# Augment data with spbp model information

Add fitted values, residuals, and optional survival predictions to the
training (or supplied) data.

## Usage

``` r
# S3 method for class 'spbp'
augment(
  x,
  data = NULL,
  eval_time = NULL,
  type = c("martingale", "deviance", "cox-snell", "coxsnell"),
  ...
)
```

## Arguments

- x:

  A fitted `"spbp"` object.

- data:

  Optional data frame; defaults to the training data when available.

- eval_time:

  Optional numeric vector of times for nested survival predictions (same
  structure as `predict(x, type = "survival")`).

- type:

  Residual type passed to
  [`residuals.spbp`](https://rvpanaro.github.io/spsurv/reference/residuals.spbp.md).

- ...:

  Not used.

## Value

A `data.frame` with original rows plus `.residual` and, when `eval_time`
is set, a list-column `.pred`.
