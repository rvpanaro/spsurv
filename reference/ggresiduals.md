# ggplot2 residual diagnostic plots for spbp models

Convenience wrapper around
[`residuals.spbp`](https://rvpanaro.github.io/spsurv/reference/residuals.spbp.md)
that returns a ggplot2 object for martingale, deviance, or Cox-Snell
residuals.

## Usage

``` r
ggresiduals(
  object,
  type = c("martingale", "deviance", "cox-snell", "coxsnell"),
  against = c("fitted", "index"),
  ...
)
```

## Arguments

- object:

  A fitted `"spbp"` object.

- type:

  Residual type: `"martingale"` (default), `"deviance"`, `"cox-snell"`,
  or `"coxsnell"`.

- against:

  What to plot on the x-axis: `"fitted"` (default, Cox-Snell cumulative
  hazard) or `"index"` (observation index).

- ...:

  Further arguments passed to
  [`residuals.spbp`](https://rvpanaro.github.io/spsurv/reference/residuals.spbp.md).

## Value

A `ggplot` object (requires ggplot2).

## Examples

``` r
if (FALSE) { # \dontrun{
library(spsurv)
library(ggplot2)
data(veteran, package = "survival")
fit <- bpph(Surv(time, status) ~ karno, data = veteran, degree = 4)
ggresiduals(fit, type = "martingale")
} # }
```
