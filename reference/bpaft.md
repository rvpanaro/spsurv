# Bernstein Polynomial Based Accelerated Failure Time Model

Fits the BPAFT model to time-to-event data.

## Usage

``` r
bpaft(formula, degree, data, approach = c("mle", "bayes"), ...)
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
#> Priors are ignored because the MLE approach is used.

summary(fit)
#> Bernstein Polynomial based Accelerated Failure Time model
#> Call:
#> spbp.default(formula = Surv(time, status) ~ karno + celltype, 
#>     degree = degree, data = veteran, approach = "mle", model = "aft")
#> 
#>   n= 137, number of events= 128 
#> 
#>                       coef exp(coef) se(coef)     z Pr(>|z|)    
#> karno              0.03479   1.03540  0.00489  7.12  1.1e-12 ***
#> celltypesmallcell -0.74476   0.47485  0.29085 -2.56   0.0104 *  
#> celltypeadeno     -0.89813   0.40733  0.28445 -3.16   0.0016 ** 
#> celltypelarge     -0.13139   0.87687  0.28396 -0.46   0.6436    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>                   exp(coef) exp(-coef) lower .95 upper .95
#> karno                 1.035      0.966     1.026     1.045
#> celltypesmallcell     0.475      2.106     0.269     0.840
#> celltypeadeno         0.407      2.455     0.233     0.711
#> celltypelarge         0.877      1.140     0.503     1.530
#> 
#> Likelihood ratio test= 65.1  on 4 df,   p=2e-13
#> Wald test            = 90.8  on 4 df,   p=<2e-16
```
