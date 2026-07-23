# Tidy summary for fitted spbp models

Create a broom/generics-style tidy data frame with one row per model
term.

## Usage

``` r
# S3 method for class 'spbp'
tidy(
  x,
  conf.int = FALSE,
  conf.level = 0.95,
  exponentiate = FALSE,
  component = c("coef", "baseline", "all"),
  ...
)
```

## Arguments

- x:

  A fitted `"spbp"` object.

- conf.int:

  Logical; include confidence/credible interval columns.

- conf.level:

  Interval coverage level.

- exponentiate:

  Logical; if `TRUE`, exponentiate the estimate and interval columns
  (regression coefficients only).

- component:

  Which parameters to return: `"coef"` (default), `"baseline"`
  (Bernstein `gamma` weights), or `"all"`.

- ...:

  Currently unused.

## Value

A `data.frame` with columns `term`, `estimate`, `std.error`, and when
available `statistic`, `p.value`, `conf.low`, `conf.high`. Columns that
are entirely `NA` are omitted.
