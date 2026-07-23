# spbp: The BP Based Survival Analysis Function

Semiparametric Survival Analysis Using Bernstein Polynomial

## Usage

``` r
spbp(formula, ...)
```

## Arguments

- formula:

  a [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html) response with
  event time, censoring status, and optional covariates.

- ...:

  Arguments passed to
  [`spbp.default`](https://rvpanaro.github.io/spsurv/reference/spbp.default.md),
  including `data`, `model`, `approach`, and `degree`. Further arguments
  in `...` are passed to
  [`rstan::optimizing`](https://mc-stan.org/rstan/reference/stanmodel-method-optimizing.html)
  (MLE) or
  [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)
  (Bayes), e.g. `iter`, `chains`, `init`.

## Value

An object of class `"spbp"`. See
[`spbp.default`](https://rvpanaro.github.io/spsurv/reference/spbp.default.md)
for the list of components (`coefficients`, `bp.param`, `degree`, etc.).

## Details

Fits Bernstein PH, PO, or AFT models to survival data via Stan (MLE or
Bayesian).

The generic dispatches to
[`spbp.default`](https://rvpanaro.github.io/spsurv/reference/spbp.default.md)
for formula objects. Convenience wrappers
[`bpph`](https://rvpanaro.github.io/spsurv/reference/bpph.md),
[`bppo`](https://rvpanaro.github.io/spsurv/reference/bppo.md), and
[`bpaft`](https://rvpanaro.github.io/spsurv/reference/bpaft.md) fix the
model family. See
[`vignette("getting-started", package = "spsurv")`](https://rvpanaro.github.io/spsurv/articles/getting-started.md)
for a tutorial,
[`vignette("model-families", package = "spsurv")`](https://rvpanaro.github.io/spsurv/articles/model-families.md)
for PH / PO / AFT comparison, and
[`vignette("bp-degree", package = "spsurv")`](https://rvpanaro.github.io/spsurv/articles/bp-degree.md)
for choosing the Bernstein polynomial degree.

## See also

[`spbp.default`](https://rvpanaro.github.io/spsurv/reference/spbp.default.md),
[`bpph`](https://rvpanaro.github.io/spsurv/reference/bpph.md),
[`bppo`](https://rvpanaro.github.io/spsurv/reference/bppo.md),
[`bpaft`](https://rvpanaro.github.io/spsurv/reference/bpaft.md),
[`bernstein`](https://rvpanaro.github.io/spsurv/reference/bernstein.md)

## Examples

``` r

library("spsurv")
data("veteran", package = "survival")
#> Warning: data set ‘veteran’ not found

fit_mle <- spbp(Surv(time, status) ~ karno + factor(celltype),
  data = veteran, model = "po"
)
summary(fit_mle)
#> Call:
#> spbp.default(formula = Surv(time, status) ~ karno + factor(celltype), 
#>     data = veteran, model = "po", approach = "mle", degree = 12)
#> 
#> Bernstein PO model: 
#> Regression coefficients:
#>                           Estimate    2.5%   97.5% Std. Error z value Pr(>|z|)
#> karno                      -0.0613 -0.0785 -0.0441     0.0088    -7.0    3e-12
#> factor(celltype)smallcell   1.2886  0.4316  2.1456     0.4373     2.9    0.003
#> factor(celltype)adeno       1.4380  0.5095  2.3664     0.4737     3.0    0.002
#> factor(celltype)large       0.1080 -0.8065  1.0224     0.4666     0.2    0.817
#>                              
#> karno                     ***
#> factor(celltype)smallcell ** 
#> factor(celltype)adeno     ** 
#> factor(celltype)large        
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Exponentiated coefficients:
#>                           Estimate 2.5% 97.5%
#> karno                         0.94 0.92   1.0
#> factor(celltype)smallcell     3.63 1.54   8.5
#> factor(celltype)adeno         4.21 1.66  10.7
#> factor(celltype)large         1.11 0.45   2.8
#> 
#> --- 
#> loglik = -708   AIC = 1448 

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
#> Chain 1: Gradient evaluation took 4.8e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.48 seconds.
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
#> Chain 1:  Elapsed Time: 0.197 seconds (Warm-up)
#> Chain 1:                0.143 seconds (Sampling)
#> Chain 1:                0.34 seconds (Total)
#> Chain 1: 
#> Warning: The largest R-hat is 1.08, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> 
#> Extracting posterior draws (this may take a moment)...

summary(fit_bayes)
#> Call:
#> spbp.default(formula = Surv(time, status) ~ karno + factor(celltype), 
#>     data = veteran, approach = "bayes", model = "po", priors = list(beta = c("normal(0,5)"), 
#>         gamma = "halfnormal(0,5)"), cores = 1, chains = 1, iter = 300, 
#>     degree = 12)
#> 
#> Bayesian Bernstein PO model: 
#> Regression coefficients:
#>                           Estimate   2.5%  97.5% Std. Error
#> karno                       -0.055 -0.071 -0.040        0.0
#> factor(celltype)smallcell    1.137  0.328  1.904        0.5
#> factor(celltype)adeno        1.215  0.341  1.955        0.4
#> factor(celltype)large        0.032 -0.932  0.842        0.5
#> 
#> Exponentiated coefficients:
#>                           Estimate 2.5% 97.5%
#> karno                         0.95 0.93   1.0
#> factor(celltype)smallcell     3.46 1.26   6.3
#> factor(celltype)adeno         3.70 1.11   6.8
#> factor(celltype)large         1.17 0.39   2.3
#> 
#> --- 
#> DIC = 1441   WAIC = -720 
```
