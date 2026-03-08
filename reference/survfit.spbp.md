# BP-based model survival curves

Compute survival curves for a fitted
[`spbp`](https://rvpanaro.github.io/spsurv/reference/spbp.md) model.

## Usage

``` r
# S3 method for class 'spbp'
survfit(
  formula,
  newdata,
  times,
  se.fit = TRUE,
  interval = 0.95,
  type = c("log", "log-log", "plain"),
  ...
)
```

## Arguments

- formula:

  An object of class `"spbp"` returned by
  [`spbp`](https://rvpanaro.github.io/spsurv/reference/spbp.md).

- newdata:

  Optional data frame used to obtain survival curves for specific
  covariate values.

- times:

  Optional numeric vector of time points at which to return estimates.

- se.fit:

  Logical; if `TRUE`, compute standard errors.

- interval:

  Confidence level for intervals (e.g. `0.95`).

- type:

  Character; confidence interval transformation. One of `"log"`,
  `"log-log"`, or `"plain"`.

- ...:

  Further arguments (currently ignored or reserved for future use).

## Value

An object of class `"survfit"`.

## See also

[`spbp`](https://rvpanaro.github.io/spsurv/reference/spbp.md),
[`survfit`](https://rdrr.io/pkg/survival/man/survfit.html).

## Examples

``` r
library(spsurv)
data(veteran, package = "survival")
#> Warning: data set ‘veteran’ not found
fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran)
#> Priors are ignored because the MLE approach is used.
survfit(fit)
#> Call: survfit.spbp(formula = fit)
#> 
#>        n events median 0.95LCL 0.95UCL
#> [1,] 137    128     73       1      NA
```
