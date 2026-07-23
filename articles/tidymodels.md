# tidymodels and tidybayes integration

**spsurv** registers [parsnip](https://parsnip.tidymodels.org/) engines
for censored regression when the package is loaded (after **parsnip**).
Load **censored** or **tidymodels** for survival metrics and workflows.

``` r

library(spsurv)
library(survival)
library(parsnip)
library(censored)
library(workflows)
data(veteran)
```

## parsnip engines

| Model | spsurv function | Engines |
|----|----|----|
| [`proportional_hazards()`](https://parsnip.tidymodels.org/reference/proportional_hazards.html) | `bpph` | `spsurv` (MLE), `spsurv_bayes` |
| [`proportional_odds()`](https://rvpanaro.github.io/spsurv/reference/proportional_odds.md) | `bppo` | `spsurv`, `spsurv_bayes` |
| [`survival_reg()`](https://parsnip.tidymodels.org/reference/survival_reg.html) | `bpaft` | `spsurv`, `spsurv_bayes` |

Use `scale = FALSE` when covariates are preprocessed with **recipes**
(the package default-scales internally when `scale = TRUE`).

``` r

spec <- proportional_hazards() |>
  set_engine("spsurv", degree = 5L, scale = FALSE, init = 0)

fit <- fit(spec, Surv(time, status) ~ karno + celltype, data = veteran)
predict(fit, veteran[1:2, ], type = "survival", eval_time = c(100, 200))
#> # A tibble: 2 × 1
#>   .pred       
#>   <I<list>>   
#> 1 <df [2 × 2]>
#> 2 <df [2 × 2]>
```

## workflows

``` r

wf <- workflow() |>
  add_formula(Surv(time, status) ~ karno + celltype) |>
  add_model(
    proportional_hazards() |>
      set_engine("spsurv", degree = 5L, scale = FALSE, init = 0)
  )

wf_fit <- fit(wf, data = veteran)
predict(wf_fit, veteran[1:3, ], type = "time")
#> # A tibble: 3 × 1
#>   .pred_time
#>        <dbl>
#> 1       138.
#> 2       181.
#> 3       138.
```

## Convenience constructor

[`bp_survival_reg()`](https://rvpanaro.github.io/spsurv/reference/bp_survival_reg.md)
maps `family = "ph"`, `"po"`, or `"aft"` to the appropriate parsnip
specification.

``` r

bp_survival_reg(family = "ph", engine = "spsurv")
#> Proportional Hazards Model Specification (censored regression)
#> 
#> Computational engine: spsurv
```

## Censored predictions from `spbp` fits

Direct fits support censored-style
[`predict()`](https://rdrr.io/r/stats/predict.html) types (same
structure as **censored**):

``` r

fit <- bpph(Surv(time, status) ~ karno, data = veteran, approach = "mle", init = 0)
predict(fit, veteran[1:2, ], type = "survival", eval_time = c(50, 100))
#> # A tibble: 2 × 1
#>   .pred       
#>   <I<list>>   
#> 1 <df [2 × 2]>
#> 2 <df [2 × 2]>
predict(fit, veteran[1:2, ], type = "linear_pred")
#> # A tibble: 2 × 1
#>   .pred_linear_pred
#>               <dbl>
#> 1             -2.03
#> 2             -2.37
generics::augment(fit, data = veteran[1:5, ])
#>   trt celltype time status karno diagtime age prior  .residual
#> 1   1 squamous   72      1    60        7  69     0  0.3107689
#> 2   1 squamous  411      1    70        5  64    10 -1.3920945
#> 3   1 squamous  228      1    60        3  38     0 -0.9255946
#> 4   1 squamous  126      1    60        9  63    10 -0.1872495
#> 5   1 squamous  118      1    70       11  65    10  0.2024279
```

Curve predictions (default) are unchanged:
`predict(fit, times = seq(0, 200, 2))`.

## Bayesian: tidybayes

For `approach = "bayes"`, use **posterior** and **tidybayes** after
fitting:

``` r

fit_bayes <- bpph(
  Surv(time, status) ~ karno,
  data = veteran,
  approach = "bayes",
  degree = 4L,
  iter = 200,
  warmup = 100,
  chains = 1,
  cores = 1,
  init = 0
)
#> 
#> SAMPLING FOR MODEL 'spbp' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 4.1e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.41 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: WARNING: There aren't enough warmup iterations to fit the
#> Chain 1:          three stages of adaptation as currently configured.
#> Chain 1:          Reducing each adaptation stage to 15%/75%/10% of
#> Chain 1:          the given number of warmup iterations:
#> Chain 1:            init_buffer = 15
#> Chain 1:            adapt_window = 75
#> Chain 1:            term_buffer = 10
#> Chain 1: 
#> Chain 1: Iteration:   1 / 200 [  0%]  (Warmup)
#> Chain 1: Iteration:  20 / 200 [ 10%]  (Warmup)
#> Chain 1: Iteration:  40 / 200 [ 20%]  (Warmup)
#> Chain 1: Iteration:  60 / 200 [ 30%]  (Warmup)
#> Chain 1: Iteration:  80 / 200 [ 40%]  (Warmup)
#> Chain 1: Iteration: 100 / 200 [ 50%]  (Warmup)
#> Chain 1: Iteration: 101 / 200 [ 50%]  (Sampling)
#> Chain 1: Iteration: 120 / 200 [ 60%]  (Sampling)
#> Chain 1: Iteration: 140 / 200 [ 70%]  (Sampling)
#> Chain 1: Iteration: 160 / 200 [ 80%]  (Sampling)
#> Chain 1: Iteration: 180 / 200 [ 90%]  (Sampling)
#> Chain 1: Iteration: 200 / 200 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.094 seconds (Warm-up)
#> Chain 1:                0.077 seconds (Sampling)
#> Chain 1:                0.171 seconds (Total)
#> Chain 1:

dr <- as_draws_df.spbp(fit_bayes)
head(dr[, c(".chain", ".iteration", ".draw", "beta[karno]")])
#> # A draws_df: 6 iterations, 1 chains, and 1 variables
#>   beta[karno]
#> 1      -0.029
#> 2      -0.032
#> 3      -0.026
#> 4      -0.039
#> 5      -0.032
#> 6      -0.027
#> # ... hidden reserved variables {'.chain', '.iteration', '.draw'}
```

With **tidybayes** loaded,
[`spread_draws()`](https://mjskay.github.io/tidybayes/reference/spread_draws.html),
[`gather_draws()`](https://mjskay.github.io/tidybayes/reference/spread_draws.html),
and
[`tidy_draws()`](https://mjskay.github.io/tidybayes/reference/tidy_draws.html)
dispatch on `spbp` objects. Draw-level survival curves:

``` r

long <- spread_surv_draws.spbp(
  fit_bayes,
  times = c(50, 100, 150),
  newdata = veteran[1, , drop = FALSE]
)
head(long)
#>     trt celltype time status karno diagtime age prior .chain .iteration .draw
#> 1     1 squamous   72      1    60        7  69     0      1          1     1
#> 1.1   1 squamous   72      1    60        7  69     0      1          1     1
#> 1.2   1 squamous   72      1    60        7  69     0      1          1     1
#> 1.3   1 squamous   72      1    60        7  69     0      1          2     2
#> 1.4   1 squamous   72      1    60        7  69     0      1          2     2
#> 1.5   1 squamous   72      1    60        7  69     0      1          2     2
#>     time      surv id
#> 1     50 0.5727629  1
#> 1.1  100 0.3467720  1
#> 1.2  150 0.2210768  1
#> 1.3   50 0.6566015  1
#> 1.4  100 0.4242396  1
#> 1.5  150 0.2723965  1
```

See the *Bayesian analysis with Stan* vignette for priors and
convergence.
