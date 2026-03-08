# spbp: The BP Based Survival Analysis Function

Semiparametric Survival Analysis Using Bernstein Polynomial

## Usage

``` r
spbp(formula, ...)
```

## Arguments

- formula:

  a Surv object with time to event, status and explanatory terms.

- ...:

  Arguments passed to \`rstan::sampling\` (e.g. iter, chains) or
  \`rstan::optimizing\`.

## Value

An object of class 'spbp'.

## Details

Fits Bernstein Polynomial based Proportional regression to survival
data.

## See also

[`spbp.default`](https://rvpanaro.github.io/spsurv/reference/spbp.default.md)

[`spbp.default`](https://rvpanaro.github.io/spsurv/reference/spbp.default.md),
[`bpph`](https://rvpanaro.github.io/spsurv/reference/bpph.md),
[`bppo`](https://rvpanaro.github.io/spsurv/reference/bppo.md),
[`bpaft`](https://rvpanaro.github.io/spsurv/reference/bpaft.md),
<https://mc-stan.org/users/documentation/>

## Examples

``` r
library("spsurv")
data("veteran", package = "survival")
#> Warning: data set ‘veteran’ not found

fit_mle <- spbp(Surv(time, status) ~ karno + factor(celltype),
  data = veteran, model = "po"
)
#> Priors are ignored because the MLE approach is used.
summary(fit_mle)
#> Bernstein Polynomial based Proportional Odds model
#> Call:
#> spbp.default(formula = Surv(time, status) ~ karno + factor(celltype), 
#>     data = veteran, model = "po", approach = "mle")
#> 
#>   n= 137, number of events= 128 
#> 
#>                               coef exp(coef) se(coef)     z Pr(>|z|)    
#> karno                     -0.06126   0.94058  0.00877 -6.98  2.9e-12 ***
#> factor(celltype)smallcell  1.28840   3.62697  0.43724  2.95   0.0032 ** 
#> factor(celltype)adeno      1.43779   4.21138  0.47369  3.04   0.0024 ** 
#> factor(celltype)large      0.10787   1.11390  0.46658  0.23   0.8172    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>                           exp(coef) exp(-coef) lower .95 upper .95
#> karno                         0.941      1.063     0.925     0.957
#> factor(celltype)smallcell     3.627      0.276     1.539     8.545
#> factor(celltype)adeno         4.211      0.237     1.664    10.657
#> factor(celltype)large         1.114      0.898     0.446     2.780
#> 
#> Likelihood ratio test= 73.3  on 4 df,   p=5e-15
#> Wald test            = 65.3  on 4 df,   p=2e-13

fit_bayes <- spbp(Surv(time, status) ~ karno + factor(celltype),
  data = veteran, model = "po", approach = "bayes",
  cores = 1, iter = 300, chains = 1,
  priors = list(
    beta = c("normal(0,5)"),
    gamma = "halfnormal(0,5)"
  )
)
#> 
#> SAMPLING FOR MODEL 'spbp' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 5.6e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.56 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 300 [  0%]  (Warmup)
#> Chain 1: Iteration:  30 / 300 [ 10%]  (Warmup)
#> Chain 1: Iteration:  60 / 300 [ 20%]  (Warmup)
#> Chain 1: Iteration:  90 / 300 [ 30%]  (Warmup)
#> Chain 1: Iteration: 120 / 300 [ 40%]  (Warmup)
#> Chain 1: Iteration: 150 / 300 [ 50%]  (Warmup)
#> Chain 1: Iteration: 151 / 300 [ 50%]  (Sampling)
#> Chain 1: Iteration: 180 / 300 [ 60%]  (Sampling)
#> Chain 1: Iteration: 210 / 300 [ 70%]  (Sampling)
#> Chain 1: Iteration: 240 / 300 [ 80%]  (Sampling)
#> Chain 1: Iteration: 270 / 300 [ 90%]  (Sampling)
#> Chain 1: Iteration: 300 / 300 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.173 seconds (Warm-up)
#> Chain 1:                0.128 seconds (Sampling)
#> Chain 1:                0.301 seconds (Total)
#> Chain 1: 
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> 
#> Extracting posterior draws (this may take a moment)...

summary(fit_bayes)
#> Bayesian Bernstein Polynomial based Proportional Odds model
#> Call:
#> spbp.default(formula = Surv(time, status) ~ karno + factor(celltype), 
#>     data = veteran, approach = "bayes", model = "po", priors = list(beta = c("normal(0,5)"), 
#>         gamma = "halfnormal(0,5)"), cores = 1, chains = 1, iter = 300)
#> 
#>   n= 137, number of events= 128 
#> 
#>                           mean(coef) median(coef) mean(exp(coef)) sd(coef)
#> karno                        -0.0564      -0.0576          0.9452     0.01
#> factor(celltype)smallcell     1.1890       1.1687          3.5626     0.41
#> factor(celltype)adeno         1.3116       1.3479          4.0430     0.41
#> factor(celltype)large         0.1011       0.0812          1.2273     0.45
#> ---
#>                           mean(exp(coef)) mean(exp(-coef)) lower .95HPD
#> karno                               0.945            1.058        0.931
#> factor(celltype)smallcell           3.563            0.332        1.037
#> factor(celltype)adeno               4.043            0.293        1.747
#> factor(celltype)large               1.227            1.001        0.304
#>                           upper .95HPD
#> karno                            0.961
#> factor(celltype)smallcell        6.170
#> factor(celltype)adeno            8.039
#> factor(celltype)large            2.321
#> 
#> DIC=  1441   WAIC=  -720  LPML= -720 
```
