## spsurv <a href='https://rvpanaro.github.io/spsurv'><img src="https://raw.githubusercontent.com/rvpanaro/spsurv/master/inst/figures/logo.png" align="right" height="139"/></a>

[![CRAN status](https://www.r-pkg.org/badges/version/spsurv)](https://cran.r-project.org/package=spsurv) [![CRAN downloads](https://cranlogs.r-pkg.org/badges/spsurv)](https://cran.r-project.org/package=spsurv)  [![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/spsurv)](https://cran.r-project.org/package=spsurv)[![R-CMD-check](https://github.com/rvpanaro/spsurv/workflows/R-CMD-check/badge.svg)](https://github.com/rvpanaro/spsurv/actions) [![Codecov test coverage](https://codecov.io/gh/rvpanaro/spsurv/branch/master/graph/badge.svg)](https://codecov.io/gh/rvpanaro/spsurv)

`spsurv` provides flexible **semi-parametric survival regression** for right-censored data: **proportional hazards (PH)**, **proportional odds (PO)**, and **accelerated failure time (AFT)** models with smooth baseline functions estimated via **Bernstein polynomials**.

### Key features

-   PH / PO / AFT regression with Bernstein polynomial baseline functions
-   Estimation via **maximum likelihood** or **Bayesian inference**
-   Bayesian fitting via **Stan** (NUTS/HMC), with multiple prior options
-   Extensible interface for user-defined modeling workflows

## Installation

### From CRAN

-   Install the CRAN version:

``` r
install.packages("spsurv")
```

### Development version

-   Installation using the devtools package:

``` r
install.packages("remotes")
remotes::install_github("rvpanaro/spsurv")
```

## Quick start

``` r
library(survival)
library(KMsurv)
library(spsurv)

data("larynx")  # time = follow-up time, delta = event indicator

fit <- spbp(Surv(time, delta) ~ age + factor(stage), model = "ph", data = larynx)
summary(fit)  
```

## Bayesian example

``` r
set.seed(1)

fit_bayes <- spbp(
  Surv(time, delta) ~ age + factor(stage),
  model = "ph",
  approach = "bayes",
  data = larynx,
  iter = 2000,
  warmup = 1000,
  chains = 4
)

summary(fit_bayes)
```

More examples: <https://rvpanaro.github.io/spsurv/reference/index.html>

## Citation

```         
citation("spsurv")
```

## Getting help / reporting bugs

-   Bug reports and feature requests: <https://github.com/rvpanaro/spsurv/issues>
-   Questions: open a GitHub Discussion (if enabled) or file an issue with a reproducible example
-   Contact: [rvpanaro\@gmail.com](mailto:rvpanaro@gmail.com){.email}
