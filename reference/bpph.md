# Bernstein Polynomial Based Proportional Hazards Model

Fits the BPPH model to time-to-event data.

## Usage

``` r
bpph(formula, degree, data, approach = c("mle", "bayes"), ...)
```

## Arguments

- formula:

  a Surv object with time to event observations, right censoring status
  and explanatory terms.

- degree:

  Bernstein polynomial degree.

- data:

  a data.frame object.

- approach:

  Bayesian or maximum likelihood estimation methods, default is approach
  = "mle".

- ...:

  further arguments passed to or from other methods

## Value

An object of class \`spbp\`.

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
#> Priors are ignored because the MLE approach is used.

summary(fit)
#> Bernstein Polynomial based Proportional Hazards model
#> Call:
#> spbp.default(formula = Surv(time, status) ~ karno + factor(celltype), 
#>     degree = degree, data = veteran, approach = "mle", model = "ph")
#> 
#>   n= 137, number of events= 128 
#> 
#>                               coef exp(coef) se(coef)     z Pr(>|z|)    
#> karno                     -0.03101   0.96947  0.00516 -6.01  1.9e-09 ***
#> factor(celltype)smallcell  0.73446   2.08436  0.25329  2.90  0.00373 ** 
#> factor(celltype)adeno      1.13852   3.12213  0.29366  3.88  0.00011 ***
#> factor(celltype)large      0.32857   1.38898  0.27621  1.19  0.23421    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>                           exp(coef) exp(-coef) lower .95 upper .95
#> karno                         0.969       1.03     0.960     0.979
#> factor(celltype)smallcell     2.084       0.48     1.269     3.424
#> factor(celltype)adeno         3.122       0.32     1.756     5.552
#> factor(celltype)large         1.389       0.72     0.808     2.387
#> 
#> Likelihood ratio test= 59.9  on 4 df,   p=3e-12
#> Wald test            = 61.9  on 4 df,   p=1e-12
```
