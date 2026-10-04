# Bernstein AFT Model

Fits the BPAFT model to time-to-event data.

## Usage

``` r
bpaft(
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

  a Surv object with time to event observations, right censoring status
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
[`bppo`](https://rvpanaro.github.io/spsurv/reference/bppo.md) for other
BP based models.

## Examples

``` r

library("spsurv")
data("veteran", package = "survival")
#> Warning: data set ‘veteran’ not found

fit <- bpaft(Surv(time, status) ~ karno + celltype,
  data = veteran
)

summary(fit)
#> Call:
#> bpaft(formula = Surv(time, status) ~ karno + celltype, data = veteran, 
#>     approach = "mle", model = "aft")
#> 
#> Bernstein AFT model: 
#> Regression coefficients:
#>                   Estimate    2.5%   97.5% Std. Error z value Pr(>|z|)    
#> karno               0.0347  0.0251  0.0443     0.0049     7.1    1e-12 ***
#> celltypesmallcell  -0.7493 -1.3149 -0.1836     0.2886    -2.6    0.009 ** 
#> celltypeadeno      -0.9003 -1.4566 -0.3440     0.2838    -3.2    0.002 ** 
#> celltypelarge      -0.1343 -0.6891  0.4206     0.2831    -0.5    0.635    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Exponentiated coefficients:
#>                   Estimate 2.5% 97.5%
#> karno                 1.04 1.03   1.0
#> celltypesmallcell     0.47 0.27   0.8
#> celltypeadeno         0.41 0.23   0.7
#> celltypelarge         0.87 0.50   1.5
#> 
#> --- 
#> loglik = -710   AIC = 1451 
```
