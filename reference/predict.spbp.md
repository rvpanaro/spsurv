# Predicted survival as a tidy data frame

Survival (and optional CI) on a time grid, as a `data.frame` for
[`ggplot2::geom_line`](https://ggplot2.tidyverse.org/reference/geom_path.html).

## Usage

``` r
# S3 method for class 'spbp'
predict(
  object,
  newdata = NULL,
  times = NULL,
  eval_time = NULL,
  type = NULL,
  conf_type = NULL,
  interval = 0.95,
  interval.type = c("hpd", "quantile"),
  monotone = NULL,
  ...
)
```

## Arguments

- object:

  Fitted `"spbp"` from
  [`spbp`](https://rvpanaro.github.io/spsurv/reference/spbp.md) /
  [`bpph`](https://rvpanaro.github.io/spsurv/reference/bpph.md) / etc.

- newdata:

  Optional `data.frame` of covariate profiles (same convention as
  [`survfit.spbp`](https://rvpanaro.github.io/spsurv/reference/survfit.spbp.md)).

- times:

  Time grid. Default: a dense sequence from 0 to the maximum observed
  time (suitable for smooth
  [`ggplot2::geom_line`](https://ggplot2.tidyverse.org/reference/geom_path.html)
  / `geom_ribbon` plots). Pass a long `seq(...)` to override resolution;
  use observed event times only if stepwise curves are intended.

- eval_time:

  Evaluation times for `type = "survival"` (required).

- type:

  Prediction type. For tidymodels/censored compatibility use
  `"survival"`, `"time"`, or `"linear_pred"`. For survival curves
  (default), use `NULL`, `"curve"`, or a confidence transformation
  (`"log"`, `"log-log"`, `"plain"`).

- conf_type:

  Confidence interval transformation for curve predictions (`"log"`,
  `"log-log"`, `"plain"`). Used when `type` is `NULL`, `"curve"`, or a
  confidence transformation name.

- interval:

  Confidence level for curve predictions; passed to
  [`survfit.spbp`](https://rvpanaro.github.io/spsurv/reference/survfit.spbp.md).

- interval.type, monotone:

  Passed to
  [`survfit.spbp`](https://rvpanaro.github.io/spsurv/reference/survfit.spbp.md)
  (Bayesian fits).

- ...:

  Passed to
  [`survfit.spbp`](https://rvpanaro.github.io/spsurv/reference/survfit.spbp.md)
  for curve predictions.

## Value

For curve predictions, same structure as
[`as.data.frame.survfitbp`](https://rvpanaro.github.io/spsurv/reference/survfit.spbp.md).
For `type = "survival"`, a tibble with list-column `.pred` (elements
contain `.eval_time` and `.pred_survival`). For `type = "time"`, a
tibble with `.pred_time`. For `type = "linear_pred"`, a tibble with
`.pred_linear_pred`.

## See also

[`survfit.spbp`](https://rvpanaro.github.io/spsurv/reference/survfit.spbp.md),
[`augment.spbp`](https://rvpanaro.github.io/spsurv/reference/augment.spbp.md)

## Examples

``` r
data(veteran, package = "survival")
#> Warning: data set ‘veteran’ not found
fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
pr <- predict(fit, times = seq(0, 400, by = 2))
#> Warning: Bernstein-polynomial (gamma) information matrix is ill-conditioned (kappa = 1.04e+12, degree = 12); delta-method survival standard errors and confidence bands may be unreliable. Try refitting with a lower Bernstein degree (argument `degree` to bpph(), bppo(), bpaft(), or spbp()).
predict(fit, veteran[1:2, ], type = "survival", eval_time = c(100, 200))
#> Warning: Bernstein-polynomial (gamma) information matrix is ill-conditioned (kappa = 1.04e+12, degree = 12); delta-method survival standard errors and confidence bands may be unreliable. Try refitting with a lower Bernstein degree (argument `degree` to bpph(), bppo(), bpaft(), or spbp()).
#> # A tibble: 2 × 1
#>   .pred       
#>   <I<list>>   
#> 1 <df [2 × 2]>
#> 2 <df [2 × 2]>
predict(fit, veteran[1:2, ], type = "time")
#> Warning: Bernstein-polynomial (gamma) information matrix is ill-conditioned (kappa = 1.04e+12, degree = 12); delta-method survival standard errors and confidence bands may be unreliable. Try refitting with a lower Bernstein degree (argument `degree` to bpph(), bppo(), bpaft(), or spbp()).
#> # A tibble: 2 × 1
#>   .pred_time
#>        <dbl>
#> 1         72
#> 2        411
if (FALSE) { # \dontrun{
  ggplot2::ggplot(pr, ggplot2::aes(time, surv)) + ggplot2::geom_line()
} # }
```
