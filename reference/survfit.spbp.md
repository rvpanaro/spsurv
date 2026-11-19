# BP-based model survival curves

Compute survival curves for a fitted
[`spbp`](https://rvpanaro.github.io/spsurv/reference/spbp.md) model.

## Usage

``` r
# S3 method for class 'spbp'
survfit(
  formula,
  newdata = NULL,
  times = NULL,
  se.fit = TRUE,
  interval = 0.95,
  type = c("log", "log-log", "plain"),
  interval.type = c("hpd", "quantile"),
  monotone = NULL,
  tidy = FALSE,
  baseline = FALSE,
  ...
)

# S3 method for class 'survfitbp'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'survfitbp'
plot(
  x,
  conf.int,
  col = 1,
  lty = 1,
  lwd = 1.5,
  xlab = "Time",
  ylab = "Survival probability",
  ylim = c(0, 1),
  mark.time = FALSE,
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

  Optional evaluation times: a
  [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html) object (legacy),
  or a non-negative **numeric** vector (unique values are used; suitable
  for smooth
  [`ggplot2::geom_line`](https://ggplot2.tidyverse.org/reference/geom_path.html)
  plots via
  [`predict.spbp`](https://rvpanaro.github.io/spsurv/reference/predict.spbp.md)
  or `as.data.frame.survfitbp`).

- se.fit:

  Logical; if `TRUE`, compute standard errors.

- interval:

  Confidence level for intervals (e.g. `0.95`).

- type:

  Character; confidence interval transformation. One of `"log"`,
  `"log-log"`, or `"plain"`.

- interval.type:

  For Bayesian fits only: `"hpd"` (default) or `"quantile"`
  (equal-tailed).

- monotone:

  Logical; for Bayesian fits, enforce non-increasing credible-band
  limits over time so ribbons plot smoothly with
  [`ggplot2::geom_line`](https://ggplot2.tidyverse.org/reference/geom_path.html).
  Defaults to `TRUE` when `approach = "bayes"`.

- tidy:

  Logical; if `TRUE`, return a `data.frame` for ggplot2. instead of a
  `"survfit"` object.

- baseline:

  Logical; if `TRUE`, return the baseline survival curve \\S_0(t)\\ at
  observed event times (or `times`) with no covariate effect. When
  `tidy = TRUE`, returns a simple `data.frame` with columns `time` and
  `surv`.

- ...:

  Further arguments passed to
  [`survfit.coxph`](https://rdrr.io/pkg/survival/man/survfit.html) for
  the reference Cox object (e.g. `conf.type` is ignored; use `type`
  instead).

- x:

  Object from `survfit.spbp`.

- row.names, optional:

  Unused; included for S3 consistency.

- conf.int:

  Logical; draw pointwise interval limits. Defaults to `TRUE` for a
  single curve when limits are present.

- col, lty, lwd, xlab, ylab, ylim:

  Graphical parameters.

- mark.time:

  Ignored (Bernstein curves are continuous).

## Value

An object of class `"survfit"` (with classes `survfitbp`, `survfitcox`,
`survfit`).

A `data.frame` with columns `id` (curve index), `time`, `surv`, `lower`,
`upper`, `cumhaz`, `std.err`.

## Functions

- `as.data.frame(survfitbp)`: Tidy survival curves for
  [`ggplot2::geom_line`](https://ggplot2.tidyverse.org/reference/geom_path.html)
  / `geom_ribbon`.

- `plot(survfitbp)`: Plot Bernstein-polynomial survival as a smooth line
  (`type = "l"`), never a Kaplan–Meier step function. Pass
  `conf.int = TRUE` for pointwise bands (also drawn with lines).

## See also

[`spbp`](https://rvpanaro.github.io/spsurv/reference/spbp.md),
[`survfit`](https://rdrr.io/pkg/survival/man/survfit.html),
[`predict.spbp`](https://rvpanaro.github.io/spsurv/reference/predict.spbp.md),
`survfit.spbp` (`as.data.frame` method).

## Examples

``` r
library(spsurv)
data(veteran, package = "survival")
#> Warning: data set ‘veteran’ not found
fit <- bpph(Surv(time, status) ~ karno + factor(celltype), data = veteran)
survfit(fit)
#> Call: survfit.spbp(formula = fit)
#> 
#>        n events median 0.95LCL 0.95UCL
#> [1,] 137    128   77.5    60.3    94.7
plot(survfit(fit, times = seq(0, max(veteran$time), length.out = 80)))


data(veteran, package = "survival")
#> Warning: data set ‘veteran’ not found
fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
sf <- survfit(fit, times = seq(0, max(veteran$time), length.out = 80))
#> Warning: Bernstein-polynomial (gamma) information matrix is ill-conditioned (kappa = 1.04e+12, degree = 12); delta-method survival standard errors and confidence bands may be unreliable. Try refitting with a lower Bernstein degree (argument `degree` to bpph(), bppo(), bpaft(), or spbp()).
#> Warning: Observed information is indefinite; negative variances are set to NA.
head(as.data.frame(sf))
#>   id     time      surv lower upper    cumhaz std.err
#> 1  1  0.00000 1.0000000    NA    NA 0.0000000      NA
#> 2  1 12.64557 0.8947161    NA    NA 0.1112489      NA
#> 3  1 25.29114 0.7912094    NA    NA 0.2341926      NA
#> 4  1 37.93671 0.6947400    NA    NA 0.3642176      NA
#> 5  1 50.58228 0.6079767    NA    NA 0.4976187      NA
#> 6  1 63.22785 0.5317991    NA    NA 0.6314896      NA
```
