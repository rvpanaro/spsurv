## spsurv <a href='https://rvpanaro.github.io/spsurv'><img src="https://raw.githubusercontent.com/rvpanaro/spsurv/master/inst/figures/logo.png" align="right" width="120" height="139" alt="spsurv package logo"/></a>

[![CRAN status](https://www.r-pkg.org/badges/version/spsurv)](https://cran.r-project.org/package=spsurv) [![CRAN downloads](https://cranlogs.r-pkg.org/badges/spsurv)](https://cran.r-project.org/package=spsurv) [![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/spsurv)](https://cran.r-project.org/package=spsurv) [![R-CMD-check](https://github.com/rvpanaro/spsurv/workflows/R-CMD-check/badge.svg)](https://github.com/rvpanaro/spsurv/actions) [![Codecov](https://codecov.io/gh/rvpanaro/spsurv/branch/master/graph/badge.svg)](https://codecov.io/gh/rvpanaro/spsurv)

## Overview

**spsurv**: An R package for semi-parametric survival analysis.

The **spsurv** package was designed to contribute with a flexible set of semi-parametric survival regression options, including **proportional hazards (PH)**, **proportional odds (PO)**, and **accelerated failure time (AFT)** models for right-censored data.

The package provides:

- Survival classes (PH, PO, AFT) extensions based on a fully likelihood-based approach for either **Bayesian** or **maximum likelihood (ML)** estimation procedures
- Smooth estimates for the unknown baseline functions based on the **Bernstein polynomial (BP)**
- Integration with **Stan** for user-defined modeling
- Six distinct prior specification options in a Bayesian analysis

Stan is an open-source platform with its own language and log-probability functions for custom likelihoods and priors. Access to Stan in R is provided via **rstan**; the package uses NUTS (No-U-Turn) sampling by default for Bayesian fits.

## Installation

### From CRAN

```r
install.packages("spsurv")
```

### Development version

```r
install.packages("devtools")
devtools::install_github("rvpanaro/spsurv")
```

## Usage

Fit a BP-based survival regression PH model. If `degree` is omitted, the default
is `ceiling(sqrt(n))` where `n` is the number of rows in `data`; the resolved
value is stored on the fit as `fit$degree`.

```r
library("KMsurv")
data("larynx")
larynx$stage <- factor(larynx$stage)

library(spsurv)

fit <- bpph(Surv(time, delta) ~ age + stage, data = larynx, approach = "mle")
summary(fit)
fit$degree   # Bernstein polynomial degree used
```

Alternatively, use the `spbp` function or set the degree with `bernstein(m)`:

```r
fit <- spbp(Surv(time, delta) ~ age + stage, model = "ph", data = larynx,
            dist = bernstein(5), approach = "mle")
summary(fit)
```

Bayesian analysis with the `approach` argument:

```r
fit2 <- spbp(Surv(time, delta) ~ age + stage,
             approach = "bayes", data = larynx,
             iter = 2000, chains = 1, warmup = 1000)
summary(fit2)
```

Nested model comparison (MLE) follows the same `anova()` patterns as
`survival::survreg`:

```r
fit0 <- bpph(Surv(time, delta) ~ 1, data = larynx, approach = "mle", degree = 5)
fit1 <- bpph(Surv(time, delta) ~ age, data = larynx, approach = "mle", degree = 5)
fit2 <- bpph(Surv(time, delta) ~ age + stage, data = larynx, approach = "mle", degree = 5)

anova(fit0, fit1)
anova(fit1, fit2)
anova(fit2)   # sequential term-wise table
```

See the [reference manual](https://rvpanaro.github.io/spsurv/reference/index.html) and
vignette `vignette("getting-started", package = "spsurv")` for more examples.

### Tidy model summaries

Coefficient and model-level summaries follow the [generics](https://generics.r-lib.org/) / broom convention:

```r
library(generics)
library(ggplot2)

fit <- bpph(Surv(time, delta) ~ age + stage, data = larynx, approach = "mle")

# One row per term (hazard ratios with 95% CI)
td <- tidy(fit, conf.int = TRUE, exponentiate = TRUE)
print(td)

# One-row model summary (n, log-likelihood, global LR test, AIC, ...)
glance(fit)

# Forest plot
ggplot(td, aes(x = estimate, y = term)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "Hazard ratio", y = NULL)
```

For publication-ready tables, pipe `tidy()` output into [gt](https://gt.rstudio.com/) or similar packages.

### Baseline inference, diagnostics, and model ranking

```r
# Baseline gamma weights with Wald table (survstan-style summary)
summary(fit, show_baseline = TRUE)
print(fit, show_baseline = TRUE)   # same via print()
print(fit, bp.param = TRUE)        # baseline-only table

# Tidy baseline weights
tidy(fit, component = "baseline", conf.int = TRUE)
confint(fit, parm = names(fit$bp.param))

# Rank nested models by AIC
f0 <- bpph(Surv(time, delta) ~ 1, data = larynx, degree = 5)
f1 <- bpph(Surv(time, delta) ~ age, data = larynx, degree = 5)
rank_models(f0, f1)                # or AIC(f0, f1)

# Baseline survival S0(t) at event times
survfit(fit, baseline = TRUE, tidy = TRUE)

# ggplot2 residual plots (requires ggplot2)
ggresiduals(fit, type = "martingale")
```

If you also use the **survstan** package, qualify spsurv helpers after
`library(survstan)` — e.g. `spsurv::estimates(fit)` and `spsurv::se(fit)`.

### tidymodels and tidybayes

Load **parsnip** (or **tidymodels**) before **spsurv** so censored-regression
engines register. Use `scale = FALSE` when preprocessing with **recipes**.

```r
library(parsnip)
library(censored)
library(workflows)
library(spsurv)

wf <- workflow() |>
  add_formula(Surv(time, status) ~ karno) |>
  add_model(
    proportional_hazards() |>
      set_engine("spsurv", degree = 5, scale = FALSE, init = 0)
  )

fit_wf <- fit(wf, data = veteran)
predict(fit_wf, veteran[1:2, ], type = "survival", eval_time = 100)
```

Bayesian fits support **tidybayes** via `as_draws_df.spbp()`, `spread_draws()`,
and `spread_surv_draws.spbp()`. See `vignette("tidymodels", package = "spsurv")`.

## Troubleshooting

### Bayesian convergence checks

For Bayesian fits, always check:

- divergent transitions (target: 0)
- split-\(\hat R\) (target: close to 1, typically < 1.01)
- effective sample sizes (avoid very low bulk/tail ESS)

If warnings appear, use this escalation order:

1. increase `adapt_delta` (for example, from 0.8 to 0.9 or 0.95)
2. increase `iter` and `warmup`
3. re-check divergences, `Rhat`, and ESS before interpretation

Higher `adapt_delta` usually reduces divergences but increases runtime,
so tuning should balance stability and compute cost.

Please report issues at [https://github.com/rvpanaro/spsurv/issues](https://github.com/rvpanaro/spsurv/issues) or contact the maintainer (see DESCRIPTION).

## Links

- [Download from CRAN](https://cloud.r-project.org/package=spsurv)
- [Browse source code](https://github.com/rvpanaro/spsurv)
- [Report a bug](https://github.com/rvpanaro/spsurv/issues)

## License

GPL-3
