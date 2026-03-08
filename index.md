## spsurv [![spsurv package logo](https://raw.githubusercontent.com/rvpanaro/spsurv/master/inst/figures/logo.png)](https://rvpanaro.github.io/spsurv)

[![CRAN
status](https://www.r-pkg.org/badges/version/spsurv)](https://cran.r-project.org/package=spsurv)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/spsurv)](https://cran.r-project.org/package=spsurv)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/grand-total/spsurv)](https://cran.r-project.org/package=spsurv)
[![R-CMD-check](https://github.com/rvpanaro/spsurv/workflows/R-CMD-check/badge.svg)](https://github.com/rvpanaro/spsurv/actions)
[![Codecov](https://codecov.io/gh/rvpanaro/spsurv/branch/master/graph/badge.svg)](https://codecov.io/gh/rvpanaro/spsurv)

## Overview

**spsurv**: An R package for semi-parametric survival analysis.

The **spsurv** package was designed to contribute with a flexible set of
semi-parametric survival regression options, including **proportional
hazards (PH)**, **proportional odds (PO)**, and **accelerated failure
time (AFT)** models for right-censored data.

The package provides:

- Survival classes (PH, PO, AFT) extensions based on a fully
  likelihood-based approach for either **Bayesian** or **maximum
  likelihood (ML)** estimation procedures
- Smooth estimates for the unknown baseline functions based on the
  **Bernstein polynomial (BP)**
- Integration with **Stan** for user-defined modeling
- Six distinct prior specification options in a Bayesian analysis

Stan is an open-source platform with its own language and
log-probability functions for custom likelihoods and priors. Access to
Stan in R is provided via **rstan**; the package uses NUTS (No-U-Turn)
sampling by default for Bayesian fits.

## Installation

### From CRAN

``` r
install.packages("spsurv")
```

### Development version

``` r
install.packages("devtools")
devtools::install_github("rvpanaro/spsurv")
```

## Usage

Fit a BP-based survival regression PH model:

``` r
library("KMsurv")
data("larynx")

library(spsurv)

fit <- bpph(Surv(time, delta) ~ age + factor(stage), model = "ph", data = larynx)
summary(fit)
```

Alternatively, use the `spbp` function:

``` r
fit <- spbp(Surv(time, delta) ~ age + factor(stage), model = "ph", data = larynx)
summary(fit)
```

Bayesian analysis with the `approach` argument:

``` r
fit2 <- spbp(Surv(time, delta) ~ age + factor(stage),
             approach = "bayes", data = larynx,
             iter = 2000, chains = 1, warmup = 1000)
summary(fit2)
```

See the [reference
manual](https://rvpanaro.github.io/spsurv/reference/index.html) for more
examples.

## Troubleshooting

Please report issues at <https://github.com/rvpanaro/spsurv/issues> or
contact the maintainer (see DESCRIPTION).

## Links

- [Download from CRAN](https://cloud.r-project.org/package=spsurv)
- [Browse source code](https://github.com/rvpanaro/spsurv)
- [Report a bug](https://github.com/rvpanaro/spsurv/issues)

## License

GPL-3
