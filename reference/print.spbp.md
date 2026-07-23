# Bernstein Polynomial Based Regression Object Print

Bernstein Polynomial Based Regression Object Print

## Usage

``` r
# S3 method for class 'spbp'
print(
  x,
  bp.param = FALSE,
  show_baseline = FALSE,
  digits = 2,
  signif.stars = getOption("show.signif.stars"),
  what = "summary",
  ...
)
```

## Arguments

- x:

  an object of class spbp.

- bp.param:

  print BP parameters only (alias for baseline-only display).

- show_baseline:

  logical; append Bernstein baseline Wald table after regression
  coefficients (MLE only). See also `print(fit, bp.param = TRUE)` for
  baseline-only output.

- digits:

  number of digits to display; defaults to `2`.

- signif.stars:

  see [`getOption`](https://rdrr.io/r/base/options.html).

- what:

  character vector; any of `"summary"`, `"tidy"`, and `"glance"` for
  broom-style coefficient and model-level tables.

- ...:

  further arguments passed to
  [`summary.spbp`](https://rvpanaro.github.io/spsurv/reference/summary.spbp.md).

## Value

none
