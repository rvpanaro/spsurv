# Bernstein Polynomial Based Proportional Odds Model

Fits the BPPO model to time-to-event data.

## Usage

``` r
bppo(formula, degree, data, approach = c("mle", "bayes"), ...)
```

## Arguments

- formula:

  a Surv object with time-to-event observations, right censoring status
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
#> Priors are ignored because the MLE approach is used.

summary(fit)
#> Bernstein Polynomial based Proportional Odds model
#> Call:
#> spbp.default(formula = Surv(time, status) ~ karno + celltype, 
#>     degree = degree, data = veteran, approach = "mle", model = "po")
#> 
#>   n= 137, number of events= 128 
#> 
#>                       coef exp(coef) se(coef)     z Pr(>|z|)    
#> karno             -0.06126   0.94058  0.00877 -6.98  2.9e-12 ***
#> celltypesmallcell  1.28841   3.62703  0.43724  2.95   0.0032 ** 
#> celltypeadeno      1.43774   4.21117  0.47368  3.04   0.0024 ** 
#> celltypelarge      0.10756   1.11356  0.46658  0.23   0.8177    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>                   exp(coef) exp(-coef) lower .95 upper .95
#> karno                 0.941      1.063     0.925     0.957
#> celltypesmallcell     3.627      0.276     1.539     8.545
#> celltypeadeno         4.211      0.237     1.664    10.656
#> celltypelarge         1.114      0.898     0.446     2.779
#> 
#> Likelihood ratio test= 73.3  on 4 df,   p=5e-15
#> Wald test            = 65.4  on 4 df,   p=2e-13
```
