# Bernstein PO Model

Fits the BPPO model to time-to-event data.

## Usage

``` r
bppo(
  formula,
  degree = NULL,
  data,
  approach = c("mle", "bayes"),
  dist = NULL,
  baseline = NULL,
  ...
)
```

## Arguments

- formula:

  a Surv object with time-to-event observations, right censoring status
  and explanatory terms.

- degree:

  Bernstein polynomial degree. If omitted, defaults to
  `ceiling(sqrt(n))` for `n = nrow(data)` (see
  [`spbp`](https://rvpanaro.github.io/spsurv/reference/spbp.md)).

- data:

  a data.frame object.

- approach:

  Bayesian or maximum likelihood estimation methods, default is approach
  = "mle".

- dist:

  optional baseline specification; use
  [`bernstein`](https://rvpanaro.github.io/spsurv/reference/bernstein.md)`(m)`
  for the Bernstein polynomial degree.

- baseline:

  optional alias for `dist`.

- ...:

  further arguments passed to or from other methods

## Value

An object of class `spbp`, including component `degree`.

## See also

[`spbp`](https://rvpanaro.github.io/spsurv/reference/spbp.md),
[`bpph`](https://rvpanaro.github.io/spsurv/reference/bpph.md) and
[`bpaft`](https://rvpanaro.github.io/spsurv/reference/bpaft.md) for
other BP based models.

## Examples

``` r

library("spsurv")
data("veteran", package = "survival")
#> Warning: data set ‘veteran’ not found

fit <- bppo(Surv(time, status) ~ karno + celltype,
  data = veteran
)

summary(fit)
#> Call:
#> bppo(formula = Surv(time, status) ~ karno + celltype, data = veteran, 
#>     approach = "mle", model = "po")
#> 
#> Bernstein PO model: 
#> Regression coefficients:
#>                   Estimate    2.5%   97.5% Std. Error z value Pr(>|z|)    
#> karno              -0.0602 -0.0770 -0.0434     0.0086    -7.0    2e-12 ***
#> celltypesmallcell   1.2790  0.4198  2.1382     0.4384     2.9    0.004 ** 
#> celltypeadeno       1.4388  0.5074  2.3702     0.4752     3.0    0.002 ** 
#> celltypelarge       0.1239 -0.7910  1.0387     0.4668     0.3    0.791    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Exponentiated coefficients:
#>                   Estimate 2.5% 97.5%
#> karno                 0.94 0.93   1.0
#> celltypesmallcell     3.59 1.52   8.5
#> celltypeadeno         4.22 1.66  10.7
#> celltypelarge         1.13 0.45   2.8
#> 
#> --- 
#> loglik = -709   AIC = 1450 
```
