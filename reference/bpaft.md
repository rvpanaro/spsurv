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
#>     approach = "mle", model = "aft", degree = 12)
#> 
#> Bernstein AFT model (degree = 12): 
#> Regression coefficients:
#>                   Estimate   2.5%  97.5% Std. Error z value Pr(>|z|)    
#> karno                0.035  0.025  0.044      0.005     7.0    2e-12 ***
#> celltypesmallcell   -0.746 -1.323 -0.168      0.295    -2.5    0.011 *  
#> celltypeadeno       -0.899 -1.456 -0.342      0.284    -3.2    0.002 ** 
#> celltypelarge       -0.132 -0.690  0.426      0.284    -0.5    0.643    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Exponentiated coefficients:
#>                   Estimate 2.5% 97.5%
#> karno                 1.04 1.03   1.0
#> celltypesmallcell     0.47 0.27   0.8
#> celltypeadeno         0.41 0.23   0.7
#> celltypelarge         0.88 0.50   1.5
#> 
#> --- 
#> loglik = -710   AIC = 1451 
```
