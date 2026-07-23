# Choosing PH, PO, and AFT models

`spsurv` implements three standard survival regression frameworks under
one likelihood and syntax:

| Family | Function | Covariate effect | Typical interpretation |
|----|----|----|----|
| PH | [`bpph()`](https://rvpanaro.github.io/spsurv/reference/bpph.md) | Multiplicative on hazard | Hazard ratio |
| PO | [`bppo()`](https://rvpanaro.github.io/spsurv/reference/bppo.md) | Multiplicative on cumulative odds | Odds ratio |
| AFT | [`bpaft()`](https://rvpanaro.github.io/spsurv/reference/bpaft.md) | Additive on log-time scale | Time ratio |

## Before you model

``` r

library(spsurv)
library(survival)
library(ggplot2)
data(veteran)
veteran$celltype <- factor(veteran$celltype)
```

``` r

km_ct <- survfit(Surv(time, status) ~ celltype, data = veteran)
km_long <- data.frame(
  time = km_ct$time,
  surv = km_ct$surv,
  celltype = rep(levels(veteran$celltype), km_ct$strata)
)
ggplot(km_long, aes(x = time, y = surv, color = celltype)) +
  geom_step(linewidth = 0.6) +
  labs(x = "Time (days)", y = "Survival probability", color = "Cell type") +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
```

![Kaplan-Meier by cell
type.](model-families_files/figure-html/eda-km-1.png)

Kaplan-Meier by cell type.

``` r

ggplot(veteran, aes(x = factor(status), y = karno, fill = factor(status))) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.alpha = 0.4) +
  scale_x_discrete(labels = c("Censored", "Event")) +
  labs(x = NULL, y = "Karnofsky score", fill = NULL) +
  theme_bw() +
  theme(legend.position = "none")
```

![Karnofsky score by event
status.](model-families_files/figure-html/eda-karno-1.png)

Karnofsky score by event status.

The paper examples sometimes restrict to `prior == 0`; we use the full
veteran dataset here for simplicity.

### Why this is insightful

KM by cell type reveals crossing and separation patterns. Karnofsky
differs between events and censoring — a strong candidate covariate.

### What to decide

- Crossing KM curves → PO or AFT may fit better than PH.
- Karnofsky separation → include in all three family fits for fair
  comparison.

## Fit and summarize

``` r

fml <- Surv(time, status) ~ karno + celltype
deg <- 5L
```

``` r

fit_ph <- bpph(fml, degree = deg, data = veteran, approach = "mle", init = 0)
```

``` r

fit_po <- bppo(fml, degree = deg, data = veteran, approach = "mle", init = 0)
```

``` r

fit_aft <- bpaft(fml, degree = 3, data = veteran, approach = "mle", init = 0)
```

AFT uses `degree = 3` here so the Bernstein `gamma` block stays
numerically stable on veteran; PH and PO use `degree = 5`. The overlay
below plots **point estimates** of survival; when
`gamma_information_stable` is `FALSE`, bands are still computed but
[`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html) warns that
they may be unreliable.

``` r

summary(fit_ph)$coefficients[, c("coef", "exp(coef)", "Pr(>|z|)")]
#>                          coef exp(coef)     Pr(>|z|)
#> karno             -0.03061037 0.9698534 1.951228e-09
#> celltypesmallcell  0.73804877 2.0918498 3.294217e-03
#> celltypeadeno      1.15357958 3.1695182 8.086335e-05
#> celltypelarge      0.33021792 1.3912713 2.317345e-01
```

## Visual check — survival overlay

Same covariate profiles under PH, PO, and AFT:

``` r

newdata <- data.frame(karno = c(30, 70), celltype = "squamous")
plot_times <- seq(0, max(veteran$time), length.out = 121)

pr_ph <- predict(fit_ph, newdata = newdata, times = plot_times)
pr_po <- predict(fit_po, newdata = newdata, times = plot_times)
pr_aft <- predict(fit_aft, newdata = newdata, times = plot_times)

pr_ph$family <- "PH"
pr_po$family <- "PO"
pr_aft$family <- "AFT"
pr_all <- rbind(pr_ph, pr_po, pr_aft)
pr_all$label <- paste0(
  pr_all$family, " (Karnofsky ",
  newdata$karno[match(as.character(pr_all$id), as.character(seq_len(nrow(newdata))))],
  ")"
)

ggplot(pr_all, aes(x = time, y = surv, color = label, linetype = family)) +
  geom_line(linewidth = 0.5) +
  labs(x = "Time (days)", y = "Survival probability", color = NULL, linetype = "Family") +
  theme_bw() +
  theme(legend.position = "bottom")
```

![Predicted survival: PH vs PO vs AFT (Karnofsky 30 and 70,
squamous).](model-families_files/figure-html/survival-overlay-1.png)

Predicted survival: PH vs PO vs AFT (Karnofsky 30 and 70, squamous).

``` r

aic_tbl <- data.frame(
  model = c("PH", "PO", "AFT"),
  AIC = c(AIC(fit_ph), AIC(fit_po), AIC(fit_aft))
)
aic_tbl
#>   model      AIC
#> 1    PH 1450.026
#> 2    PO 1438.685
#> 3   AFT 1464.954
```

### Why this is insightful

The same covariates on different regression scales produce **different
survival curves**. AIC summarizes fit but does not replace clinical
judgment when curves diverge meaningfully.

### What to decide

- PH vs PO curves similar → PH may suffice for reporting hazard ratios.
- Large divergence or known crossing in KM → prefer PO or check
  diagnostics.
- AFT when **time ratios** are the estimand of interest.
- Lower AIC ≠ automatic winner if curves disagree clinically.

## When to use which family

- **PH** — constant hazard ratios over time.
- **PO** — proportionality on cumulative odds; useful when hazard ratios
  are not stable.
- **AFT** — direct interpretation on the time scale.

All three share the same Bernstein baseline machinery and Stan backend.
