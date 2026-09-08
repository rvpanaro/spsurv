# Bernstein PH Model

Fits the BPPH model to time-to-event data.

## Usage

``` r
bpph(
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
[`bppo`](https://rvpanaro.github.io/spsurv/reference/bppo.md) and
[`bpaft`](https://rvpanaro.github.io/spsurv/reference/bpaft.md) for
other BP based models.

## Examples

``` r

library("spsurv")
data("veteran", package = "survival")
#> Warning: data set ‘veteran’ not found

fit <- bpph(Surv(time, status) ~ karno + factor(celltype),
  data = veteran
)

summary(fit)
#> Call:
#> bpph(formula = Surv(time, status) ~ karno + factor(celltype), 
#>     data = veteran, approach = "mle", model = "ph")
#> 
#> Bernstein PH model: 
#> Regression coefficients:
#>                           Estimate    2.5%   97.5% Std. Error z value Pr(>|z|)
#> karno                      -0.0310 -0.0411 -0.0209     0.0052    -6.0    2e-09
#> factor(celltype)smallcell   0.7349  0.2392  1.2306     0.2529     2.9    0.004
#> factor(celltype)adeno       1.1387  0.5635  1.7139     0.2935     3.9    1e-04
#> factor(celltype)large       0.3289 -0.2122  0.8700     0.2761     1.2    0.234
#>                              
#> karno                     ***
#> factor(celltype)smallcell ** 
#> factor(celltype)adeno     ***
#> factor(celltype)large        
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Exponentiated coefficients:
#>                           Estimate 2.5% 97.5%
#> karno                         0.97 0.96   1.0
#> factor(celltype)smallcell     2.09 1.27   3.4
#> factor(celltype)adeno         3.12 1.76   5.6
#> factor(celltype)large         1.39 0.81   2.4
#> 
#> --- 
#> loglik = -714   AIC = 1460 
```
