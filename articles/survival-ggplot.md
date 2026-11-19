# Survival prediction and ggplot workflows

Unlike stepwise Cox `survfit` curves, `spsurv` evaluates survival on a
**dense time grid** using the fitted Bernstein baseline. Kaplan–Meier
stays a **step** function; Bernstein-polynomial curves are always
**smooth lines** (`geom_line`, or `plot(survfit(fit))`). This vignette
connects observed KM to model-based predictions.

## Before you model

``` r

library(spsurv)
library(KMsurv)
library(survival)
library(ggplot2)
data(larynx)
larynx$stage <- factor(larynx$stage)
```

``` r

km_stage <- survfit(Surv(time, delta) ~ stage, data = larynx)
km_long <- data.frame(
  time = km_stage$time,
  surv = km_stage$surv,
  stage = rep(levels(larynx$stage), km_stage$strata)
)
ggplot(km_long, aes(x = time, y = surv, color = stage)) +
  geom_step(linewidth = 0.6) +
  geom_point(
    data = subset(larynx, delta == 0),
    aes(x = time, y = 0, color = stage),
    shape = 3,
    size = 1.5,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  labs(x = "Time (years)", y = "Survival probability", color = "Stage") +
  theme_bw() +
  theme(legend.position = "bottom")
```

![Kaplan-Meier by stage with censoring marks
(+).](survival-ggplot_files/figure-html/eda-km-1.png)

Kaplan-Meier by stage with censoring marks (+).

We compare predictions at **age = 70** across all stages — a clinically
interpretable profile that holds age fixed while varying disease stage.

### Why this is insightful

KM curves show **unadjusted** stage separation. Censoring marks (`+`)
remind you that late follow-up estimates rest on fewer subjects.

### What to decide

- Clear stage ordering in KM → including `stage` in the model is
  justified.
- Crossing KM curves → proportional hazards may be questionable; see the
  *Choosing PH, PO, and AFT models* vignette.

## Fit and summarize

``` r

fit <- bpph(
  Surv(time, delta) ~ age + stage,
  degree = 5,
  data = larynx,
  approach = "mle",
  init = 0
)
summary(fit)
#> Call:
#> bpph(formula = Surv(time, delta) ~ age + stage, degree = 5L, 
#>     data = larynx, approach = "mle", init = 0, model = "ph")
#> 
#> Bernstein PH model (degree = 5): 
#> Regression coefficients:
#>        Estimate    2.5%   97.5% Std. Error z value Pr(>|z|)    
#> age      0.0197 -0.0085  0.0478     0.0143     1.4     0.17    
#> stage2   0.1729 -0.7325  1.0783     0.4620     0.4     0.71    
#> stage3   0.6521 -0.0450  1.3493     0.3557     1.8     0.07 .  
#> stage4   1.7777  0.9470  2.6085     0.4239     4.2    3e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Exponentiated coefficients:
#>        Estimate 2.5% 97.5%
#> age        1.02 0.99   1.0
#> stage2     1.19 0.48   2.9
#> stage3     1.92 0.96   3.9
#> stage4     5.92 2.58  13.6
#> 
#> --- 
#> loglik = -141   AIC = 299
```

## Prediction grid

``` r

plot_times <- seq(0, max(larynx$time), length.out = 121)
newdata <- data.frame(age = 70, stage = levels(larynx$stage))
newdata
#>   age stage
#> 1  70     1
#> 2  70     2
#> 3  70     3
#> 4  70     4
```

## `predict()` — tidy data frame

``` r

pr <- predict(fit, newdata = newdata, times = plot_times)
head(pr)
#>     age stage id       time      surv     lower     upper      cumhaz
#> 1    70     1  1 0.00000000 1.0000000 1.0000000 1.0000000 0.000000000
#> 1.1  70     1  1 0.08916667 0.9922944 0.9841838 1.0000000 0.007735405
#> 1.2  70     1  1 0.17833333 0.9847546 0.9695190 1.0000000 0.015362761
#> 1.3  70     1  1 0.26750000 0.9773718 0.9558600 0.9993677 0.022888184
#> 1.4  70     1  1 0.35666667 0.9701371 0.9430748 0.9979761 0.030317848
#> 1.5  70     1  1 0.44583333 0.9630423 0.9310441 0.9961401 0.037657975
#>         std.err
#> 1   0.000000000
#> 1.1 0.004187411
#> 1.2 0.007955462
#> 1.3 0.011355137
#> 1.4 0.014434892
#> 1.5 0.017240425
```

## Visual check — KM vs model overlay

**Panel A:** observed KM (**steps**) vs model prediction (**smooth**
dashed) per stage.

``` r

pr$stage <- newdata$stage[match(as.character(pr$id), as.character(seq_len(nrow(newdata))))]

ggplot() +
  geom_step(
    data = km_long,
    aes(x = time, y = surv, color = stage),
    linewidth = 0.5,
    alpha = 0.8
  ) +
  geom_line(
    data = pr,
    aes(x = time, y = surv, color = stage),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  labs(x = "Time (years)", y = "Survival probability", color = "Stage") +
  theme_bw() +
  theme(legend.position = "bottom")
```

![Observed KM (steps) vs smooth Bernstein prediction at age 70
(dashed).](survival-ggplot_files/figure-html/km-model-overlay-1.png)

Observed KM (steps) vs smooth Bernstein prediction at age 70 (dashed).

**Panel B:** model-based ribbons (adjusted prediction with uncertainty),
still drawn with `geom_line` rather than `geom_step`.

``` r

ggplot(pr, aes(x = time, y = surv, color = stage, fill = stage)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.5) +
  labs(x = "Time (years)", y = "Survival probability", color = NULL, fill = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")
```

![Smooth stage-specific Bernstein survival with 95%
intervals.](survival-ggplot_files/figure-html/ggplot-1.png)

Smooth stage-specific Bernstein survival with 95% intervals.

### Why this is insightful

Panel A connects the **smooth BP fit** to observed risk sets. Panel B
shows **covariate-adjusted** uncertainty — ribbons often widen late in
follow-up when data are sparse.

### What to decide

- Model curves systematically below KM → check covariate adjustment or
  misspecification.
- Wide bands late in follow-up → avoid over-interpreting tail survival.
- For publication → export
  [`predict()`](https://rdrr.io/r/stats/predict.html) output directly to
  your ggplot pipeline.

## `survfit()` with `tidy = TRUE`

``` r

sf <- survfit(fit, newdata = newdata, times = plot_times, tidy = TRUE)
head(sf)
#>     age stage id       time      surv     lower     upper      cumhaz
#> 1    70     1  1 0.00000000 1.0000000 1.0000000 1.0000000 0.000000000
#> 1.1  70     1  1 0.08916667 0.9922944 0.9841838 1.0000000 0.007735405
#> 1.2  70     1  1 0.17833333 0.9847546 0.9695190 1.0000000 0.015362761
#> 1.3  70     1  1 0.26750000 0.9773718 0.9558600 0.9993677 0.022888184
#> 1.4  70     1  1 0.35666667 0.9701371 0.9430748 0.9979761 0.030317848
#> 1.5  70     1  1 0.44583333 0.9630423 0.9310441 0.9961401 0.037657975
#>         std.err
#> 1   0.000000000
#> 1.1 0.004187411
#> 1.2 0.007955462
#> 1.3 0.011355137
#> 1.4 0.014434892
#> 1.5 0.017240425
```

## Bayesian fits

For `approach = "bayes"`, credible bands use posterior draws with
`interval.type = "hpd"` and `monotone = TRUE` (default for Bayes).

``` r

fit_bayes <- bpph(
  Surv(time, delta) ~ age + stage,
  degree = 5,
  data = larynx,
  approach = "bayes",
  iter = 400,
  chains = 1,
  cores = 1
)
predict(fit_bayes, newdata = newdata, times = plot_times, interval.type = "hpd")
```

## Interval types (MLE)

Log–log intervals are undefined at `S = 1` (`time = 0`), so we evaluate
at `time > 0` only:

``` r

loglog_times <- plot_times[plot_times > 0]
head(predict(fit, newdata = newdata, times = loglog_times, type = "log-log"))
#>     age stage id       time      surv     lower     upper      cumhaz
#> 1    70     1  1 0.08916667 0.9922944 0.9778986 0.9973263 0.007735405
#> 1.1  70     1  1 0.17833333 0.9847546 0.9584966 0.9944477 0.015362761
#> 1.2  70     1  1 0.26750000 0.9773718 0.9412743 0.9913813 0.022888184
#> 1.3  70     1  1 0.35666667 0.9701371 0.9258109 0.9881467 0.030317848
#> 1.4  70     1  1 0.44583333 0.9630423 0.9117639 0.9847653 0.037657975
#> 1.5  70     1  1 0.53500000 0.9560789 0.8988532 0.9812597 0.044914829
#>         std.err
#> 1   0.004187411
#> 1.1 0.007955462
#> 1.2 0.011355137
#> 1.3 0.014434892
#> 1.4 0.017240425
#> 1.5 0.019814418
```
